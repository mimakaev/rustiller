#!/usr/bin/env python3
"""rustiller: distiller project yaml -> pairs (Parquet) -> coolers, with bwa-mem3 and rooler.

    SRA run --download--> .sra --fastq-dump (streamed) | bwa-mem3 | rooler pairs parse--> lane pairs
    lanes of a library --rooler pairs merge --dedup--> library pairs + stats
    library pairs --rooler pairs view [--min-mapq] | rooler cload--> base cooler --zoomify--> mcool
    libraries of a group --rooler merge--> group cooler --zoomify--> mcool

Every task starts the moment its inputs exist and resources allow: a library is deduplicated and
binned while later lanes are still mapping. A task whose output exists is skipped, so re-running
the same command resumes an interrupted project. No FASTQ is ever written to disk.

    rustiller.py run project.yml [-o OUTDIR] [options]
    nextflow run main.nf --project project.yml --outdir OUTDIR      (the same on a cluster)
    rustiller.py sra-stream RUN.sra            (interleaved FASTQ on stdout; used by the mapper)
"""
import argparse, os, queue, re, shlex, shutil, subprocess, sys, threading, time

SRA_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}"
CIS_KB = [1, 2, 4, 10, 20, 40]


# ---------------------------------------------------------------- SRA -> interleaved FASTQ stream
def sra_spots(sra, vdb_dump):
    out = subprocess.run(shlex.split(vdb_dump) + ["-r", sra], capture_output=True, text=True, check=True).stdout
    m = re.search(r"first-row = ([\d,]+), row-count = ([\d,]+)", out)
    first, n = (int(x.replace(",", "")) for x in m.groups())
    return first, first + n - 1


def sra_stream(sra, procs, start, end, fastq_dump, vdb_dump, out, block=8 << 20):
    """Run `procs` fastq-dump processes over consecutive spot ranges and interleave their output.

    Blocks of a fixed byte size are taken from the processes in strict rotation and cut so that
    the two reads of a spot stay adjacent: the stream is the same on every run, and it is what
    `bwa mem -p` expects. Nothing touches the disk.
    """
    first, last = sra_spots(sra, vdb_dump)
    first, last = max(first, start or first), min(last, end or last)
    n = last - first + 1
    procs = max(1, min(procs, n // 1_000_000 or 1))
    bounds = [first + n * i // procs for i in range(procs)] + [last + 1]
    ps, qs = [], []

    def name_at(buf, pos):                       # name line of the record starting at pos
        return buf[pos:buf.find(b"\n", pos)]

    def reader(p, q):
        carry = b""
        while True:
            chunk = p.stdout.read(block)
            buf = carry + chunk
            if not chunk:
                if buf: q.put(buf)
                break
            lines = buf.count(b"\n")
            drop = lines % 4 + 1                 # newlines after the last whole record, plus its own
            cut = len(buf)
            for _ in range(drop if lines >= 4 else 0):
                cut = buf.rfind(b"\n", 0, cut)
            cut = cut + 1 if lines >= 4 else 0
            if cut:                              # do not separate the two reads of a spot: end the
                starts = [cut]                   # block after a second read, or before a first one
                for _ in range(2):
                    prev = starts[-1] - 1
                    for _ in range(4):
                        prev = buf.rfind(b"\n", 0, prev)
                    starts.append(prev + 1)
                    if prev < 0: break
                if len(starts) < 3 or name_at(buf, starts[1]) != name_at(buf, starts[2]) or starts[1] == starts[2]:
                    cut = starts[1]
            if cut: q.put(buf[:cut])
            carry = buf[cut:]
        q.put(None)

    for i in range(procs):
        cmd = shlex.split(fastq_dump) + ["--split-spot", "--skip-technical", "--defline-seq", "@$ac.$si",
                                         "--defline-qual", "+", "-N", str(bounds[i]), "-X", str(bounds[i + 1] - 1), "-Z", sra]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=block)
        q = queue.Queue(maxsize=3)
        threading.Thread(target=reader, args=(p, q), daemon=True).start()
        ps.append(p); qs.append(q)
    live = list(range(procs))
    while live:
        for i in list(live):
            b = qs[i].get()
            if b is None: live.remove(i)
            else: out.write(b)
    out.flush()
    bad = [p.wait() for p in ps]
    if any(bad): sys.exit(f"sra-stream: fastq-dump failed on {sra} (exit codes {bad})")


# ---------------------------------------------------------------- stats
def sum_stats(inputs, out):
    """Sum pairtools-style stats files; the summary/ fractions are recomputed from the sums."""
    order, s = [], {}
    for f in sorted(inputs):                     # the order of keys must not depend on the caller
        for line in open(f):
            k, v = line.rstrip("\n").split("\t")
            if k.startswith("summary/"): continue
            if k not in s: order.append(k); s[k] = 0
            s[k] += int(v)
    nodups = max(s.get("total_nodups", 0), 1)
    summary = [("summary/frac_cis", s.get("cis", 0) / nodups)]
    summary += [(f"summary/frac_cis_{kb}kb+", s.get(f"cis_{kb}kb+", 0) / nodups) for kb in CIS_KB if f"cis_{kb}kb+" in s]
    summary += [("summary/frac_dups", s.get("total_dups", 0) / max(s.get("total_mapped", 0), 1))]
    with open(out + ".tmp", "w") as o:
        for k in order:
            if k.startswith("chrom_freq/") and summary:
                for sk, sv in summary: o.write(f"{sk}\t{sv}\n")
                summary = []
            o.write(f"{k}\t{s[k]}\n")
        for sk, sv in summary: o.write(f"{sk}\t{sv}\n")
    os.replace(out + ".tmp", out)


# ---------------------------------------------------------------- the pipeline
class Task:
    def __init__(self, name, cmd, out, deps=(), cost=None, hold=None, release=None, prio=0):
        self.name, self.cmd, self.out, self.deps = name, cmd, out, list(deps)
        self.cost = cost or {}                   # taken at start, returned at the end
        self.hold = hold or {}                   # taken at start, kept after success ...
        self.release = release or {}             # ... until a later task returns it
        self.prio, self.state, self.seconds = prio, "pending", None


def mem_total_gb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemTotal"): return int(line.split()[1]) / 1e6
    return 16.0


def plan(a):
    """Read the project file into the plain settings both executors work from (the local one
    below, and main.nf through `rustiller.py plan`)."""
    import yaml
    cfg = yaml.safe_load(open(a.project))
    inp, genome = cfg["input"], cfg["input"]["genome"]
    extra = cfg.get("rustiller") or {}
    tool = lambda k: os.environ.get("RUSTILLER_" + k.upper().replace("-", "_")) or (extra.get("tools") or {}).get(k) or k
    P = {"assembly": genome["assembly_name"], "chromsizes": genome["chrom_sizes_path"],
         "tools": {k: tool(k) for k in ["rooler", "bwa-mem3", "fastq-dump", "vdb-dump"]}}

    P["index"] = index = a.index or genome.get("bwa_mem3_index_path") or re.sub(r"\.?\*$", "", genome.get("bwa_index_wildcard_path", ""))
    if not os.path.exists(index + ".bwt.2bit.64"):
        sys.exit(f"no bwa-mem3 index at {index}.* (a classic bwa index does not work). Build one with "
                 f"`bwa-mem3 index <genome.fa>` and pass --index or set input.genome.bwa_mem3_index_path")

    # map / parse / dedup options
    mp = cfg.get("map") or {}
    if (mp.get("trim_options") or "").strip(): sys.exit("map.trim_options is not supported")
    P["map_opts"] = mp.get("mapping_options") or ""
    popts = (cfg.get("parse") or {}).get("parsing_options") or ""
    if not re.search(r"--walks-policy\s+all\b", popts):
        print(f"WARNING: parse.parsing_options is {popts!r}; rooler implements `--walks-policy all` only "
              f"(every ligation junction of a read pair is reported)", file=sys.stderr)
    P["parse_opts"] = " ".join(m.group(0) for m in re.finditer(r"--(min-mapq|max-inter-align-gap|max-molecule-size)[ =]\d+", popts)).replace("=", " ")
    P["max_mismatch"] = (cfg.get("dedup") or {}).get("max_mismatch_bp", 0)

    # binning
    b = cfg.get("bin") or {}
    res = a.resolutions or b.get("resolutions") or [1000]
    if isinstance(res, str):
        if res.startswith("ladder:"):
            base = int(res.split(":")[1]); res = [base << k for k in range(40) if base << k <= 10_000_000]
        else: res = [int(x) for x in res.split(",")]
    P["resolutions"] = res = sorted(set(res))
    if any(r % res[0] for r in res): sys.exit(f"every resolution must be a multiple of the finest one ({res[0]})")
    P["balance"] = bool(b.get("balance", True))
    if (b.get("balance_options") or "").strip(): print("WARNING: bin.balance_options is ignored (rooler defaults are used)", file=sys.stderr)
    filters = {}
    for name, expr in (b.get("filters") or {"no_filter": ""}).items():
        m = re.fullmatch(r"\(?\s*mapq1\s*>=\s*(\d+)\s*\)?\s+and\s+\(?\s*mapq2\s*>=\s*(\d+)\s*\)?", (expr or "").strip())
        if not (expr or "").strip(): filters[name] = 0
        elif m and m.group(1) == m.group(2): filters[name] = int(m.group(1))
        else: sys.exit(f"bin.filters.{name}: only '(mapq1>=N) and (mapq2>=N)' and '' are supported, got {expr!r}")
    if a.filters: filters = {k: v for k, v in filters.items() if k in a.filters.split(",")}
    P["filters"] = filters

    # lanes: an SRA run (optionally a spot range of it) or two local FASTQ files
    P["lanes"] = []
    for lib, lanes in inp["raw_reads_paths"].items():
        for lane, src in lanes.items():
            m = re.fullmatch(r"sra:(\w+)(?:\?start=(\d+)&end=(\d+))?", src[0]) if len(src) == 1 else None
            if m:
                P["lanes"].append({"library": lib, "lane": lane, "sra": m.group(1),
                                   "start": int(m.group(2) or 0), "end": int(m.group(3) or 0)})
            elif len(src) == 2 and not any("://" in x or x.startswith("sra:") for x in src):
                P["lanes"].append({"library": lib, "lane": lane, "fastqs": [os.path.abspath(x) for x in src]})
            else:
                sys.exit(f"{lib}/{lane}: expected one 'sra:<run>' entry or two local FASTQ paths, got {src}")
    P["groups"] = {g: list(m) for g, m in (inp.get("library_groups") or {}).items()}
    P["dirs"] = {k: ((cfg.get("output") or {}).get("dirs") or {}).get(k, "results/" + k).rstrip("/") for k in
                 ["mapped_parsed_sorted_chunks", "pairs_library", "coolers_library", "coolers_library_group", "stats_library_group"]}
    P["dirs"]["stats_library"] = "results/stats_library"
    return P


def build(a):
    P = plan(a)
    asm, chromsizes, index, res, filters = P["assembly"], P["chromsizes"], P["index"], P["resolutions"], P["filters"]
    ROOLER, BWA, FQD, VDB = (P["tools"][k] for k in ["rooler", "bwa-mem3", "fastq-dump", "vdb-dump"])
    map_opts, parse_opts, mm, base = P["map_opts"], P["parse_opts"], P["max_mismatch"], res[0]
    balance = "--balance" if P["balance"] else ""
    me = f"{shlex.quote(sys.executable)} {shlex.quote(os.path.abspath(__file__))}"

    # directories
    OUT = os.path.abspath(a.outdir)
    TMP = os.path.abspath(a.tmpdir or os.path.join(OUT, "tmp"))
    D = {k: os.path.join(OUT, v) for k, v in P["dirs"].items()}
    LOGS, SRA = os.path.join(OUT, "logs"), os.path.join(TMP, "sra")

    # resources: one mapper (bwa-mem3 scales to all cores in one process; several smaller ones
    # are no faster) plus helpers on the remaining and oversubscribed threads
    cpus = a.cpus or len(os.sched_getaffinity(0))
    bwa_t = max(1, cpus * 3 // 4)
    parse_t = max(1, min(4, cpus // 8))
    help_t = max(1, min(8, cpus // 4))
    mem = a.mem or mem_total_gb() * 0.8
    help_mem = max(2, min(16, mem / 6))
    resources = {"cpu": cpus + help_t, "mem": mem, "net": a.downloads, "sra": a.sra_ahead}
    MAP = {"cpu": cpus, "mem": min(24, mem * 0.5)}
    HELP = {"cpu": help_t, "mem": help_mem}
    hopt = f"--threads {help_t} --mem {help_mem:.0f}"

    tasks = {}
    def add(name, *args, **kw):
        tasks[name] = Task(name, *args, **kw); return name

    def coolers(name, pairs_or_cools, deps, d, from_pairs):
        """base cooler (from pairs, or merged from library coolers) -> mcool, for every filter"""
        for fname, q in filters.items():
            cool = f"{d}/{name}.{asm}.{fname}.{base}.cool"
            if from_pairs:
                cmd = (f"{ROOLER} pairs view {pairs_or_cools} --min-mapq {q} --chroms {chromsizes} "
                       f"| {ROOLER} cload - {base} {cool}.tmp --assembly {asm} {hopt} && mv {cool}.tmp {cool}")
                dd = deps
            else:
                ins = " ".join(f"{D['coolers_library']}/{l}.{asm}.{fname}.{base}.cool" for l in pairs_or_cools)
                cmd = f"{ROOLER} merge {cool}.tmp {ins} --assembly {asm} {hopt} && mv {cool}.tmp {cool}"
                dd = [f"cload.{l}.{fname}" for l in pairs_or_cools]
            c = add(f"{'cload' if from_pairs else 'merge'}.{name}.{fname}", cmd, cool, dd, HELP, prio=4 if from_pairs else 5)
            mcool = cool[:-5] + ".mcool"
            add(f"zoomify.{name}.{fname}",
                f"{ROOLER} zoomify {cool} {mcool}.tmp --resolutions {','.join(map(str, res))} {balance} "
                f"--assembly {asm} {hopt} && mv {mcool}.tmp {mcool}", mcool, [c], HELP, prio=4 if from_pairs else 5)

    uses, by_lib = {}, {}                        # a run split over several lanes is kept until the end
    for L in P["lanes"]:
        by_lib.setdefault(L["library"], []).append(L)
        if "sra" in L: uses[L["sra"]] = uses.get(L["sra"], 0) + 1
    for lib, lanes in by_lib.items():
        lane_tasks, lane_files = [], []
        for L in lanes:
            lane = L["lane"]
            pq = f"{D['mapped_parsed_sorted_chunks']}/{lib}.{lane}.{asm}.pairs.parquet"
            part = f"{TMP}/{lib}.{lane}.parquet"
            deps, release, clean = [], {}, ""
            if "sra" in L:
                acc, rng = L["sra"], (f" --start {L['start']} --end {L['end']}" if L["end"] else "")
                sra = f"{SRA}/{acc}.sra"
                if not os.path.exists(pq):
                    deps = [f"fetch.{acc}"] if f"fetch.{acc}" in tasks else [add(f"fetch.{acc}",
                                f"curl -sS -f --retry 8 --retry-delay 10 -C - -o {sra}.part {SRA_URL.format(acc=acc)} && mv {sra}.part {sra}",
                                sra, cost={"net": 1}, hold={"sra": 1}, prio=1)]
                    release = {"sra": 1} if uses[acc] == 1 else {}
                reads = f"{me} sra-stream {sra} --procs {a.dump_procs} --fastq-dump {shlex.quote(FQD)} --vdb-dump {shlex.quote(VDB)}{rng} |"
                bwa_in = "-p -"
                if not a.keep_sra and uses[acc] == 1: clean = f" && rm -f {sra}"
            else:
                reads, bwa_in = "", " ".join(shlex.quote(x) for x in L["fastqs"])
            m = add(f"map.{lib}.{lane}",
                    f"{reads} {BWA} mem -SP5M {map_opts} -t {bwa_t} {index} {bwa_in} "
                    f"| {ROOLER} pairs parse - {part} --chroms-path {chromsizes} --assembly {asm} {parse_opts} "
                    f"--threads {parse_t} --mem 8 && mv {part} {pq}{clean}",
                    pq, deps, MAP, release=release, prio=2)
            lane_tasks.append(m); lane_files.append(pq)
        libpq = f"{D['pairs_library']}/{lib}.{asm}.nodups.pairs.parquet"
        dd = add(f"dedup.{lib}",
                 f"{ROOLER} pairs merge {libpq}.tmp {' '.join(lane_files)} --dedup {mm} "
                 f"--stats {D['stats_library']}/{lib}.{asm}.dedup.stats {hopt} && mv {libpq}.tmp {libpq}",
                 libpq, lane_tasks, HELP, prio=3)
        coolers(lib, libpq, [dd], D["coolers_library"], True)

    for group, members in P["groups"].items():
        coolers(group, members, None, D["coolers_library_group"], False)
        gstats = f"{D['stats_library_group']}/{group}.{asm}.stats"
        stats = " ".join(f"{D['stats_library']}/{l}.{asm}.dedup.stats" for l in members)
        add(f"stats.{group}", f"{me} sum-stats {gstats} {stats}", gstats, [f"dedup.{l}" for l in members], {"cpu": 1}, prio=5)

    return tasks, resources, list(D.values()) + [LOGS, TMP, SRA], LOGS


def run_pipeline(a):
    tasks, free, dirs, LOGS = build(a)
    for t in tasks.values():
        if os.path.exists(t.out): t.state = "done"
    if a.dry_run:
        for t in tasks.values():
            print(f"[{t.state:7}] {t.name}  <- {', '.join(t.deps) or '-'}\n    {t.cmd}")
        print(f"{len(tasks)} tasks; resources {free}")
        return 0
    for d in dirs: os.makedirs(d, exist_ok=True)
    lock, t_start = threading.Condition(), time.time()
    def log(msg):
        print(f"[{time.strftime('%H:%M:%S')} +{(time.time() - t_start) / 60:6.1f}m] {msg}", flush=True)

    def run(t):
        t0 = time.time()
        with open(os.path.join(LOGS, t.name + ".log"), "w") as lf:
            lf.write(t.cmd + "\n"); lf.flush()
            rc = subprocess.call(["bash", "-o", "pipefail", "-c", t.cmd], stdout=lf, stderr=lf)
        with lock:
            ok = rc == 0 and os.path.exists(t.out)
            t.state, t.seconds = ("done" if ok else "failed"), time.time() - t0
            for r, n in t.cost.items(): free[r] += n
            for r, n in t.hold.items():
                if not ok: free[r] += n
            for r, n in t.release.items(): free[r] += n
            log(f"{'done  ' if ok else 'FAILED'} {t.name}  ({t.seconds:.0f}s)")
            lock.notify_all()

    with lock:
        for t in tasks.values():                 # a downloaded run awaiting its mapper holds a slot
            if t.state == "done":
                for r, n in t.hold.items(): free[r] -= n
        log(f"{sum(t.state == 'done' for t in tasks.values())} of {len(tasks)} tasks already done; resources {free}")
        while True:
            for t in tasks.values():
                if t.state == "pending" and any(tasks[d].state in ("failed", "skipped") for d in t.deps):
                    t.state = "skipped"; log(f"skipped {t.name} (a dependency failed)")
            ready = [t for t in tasks.values() if t.state == "pending" and all(tasks[d].state == "done" for d in t.deps)]
            maps = [t for t in tasks.values() if t.name.startswith("map.")]
            waiting = [t for t in maps if t.state == "pending"]
            reserve = waiting[0].cost if waiting and not any(t.state == "running" for t in maps) else {}
            for t in sorted(ready, key=lambda t: -t.prio):   # stable: project order within a priority
                need = {**t.cost, **t.hold}
                # the mapper is the bottleneck: helpers never take the threads it is waiting for
                room = reserve if not t.name.startswith("map.") else {}
                if all(free[r] - room.get(r, 0) >= n for r, n in need.items()):
                    for r, n in need.items(): free[r] -= n
                    t.state = "running"; log(f"start  {t.name}")
                    threading.Thread(target=run, args=(t,), daemon=True).start()
            if not any(t.state in ("pending", "running") for t in tasks.values()): break
            lock.wait()

    bad = [t.name for t in tasks.values() if t.state != "done"]
    kinds = {}
    for t in tasks.values():
        if t.seconds is not None:
            k = t.name.split(".")[0]; kinds[k] = kinds.get(k, 0) + t.seconds
    log("task seconds by kind: " + ", ".join(f"{k} {v:.0f}" for k, v in kinds.items()))
    log(f"finished: {len(tasks) - len(bad)} done, {len(bad)} failed/skipped {bad}")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run a distiller project file")
    r.add_argument("project")
    r.add_argument("-o", "--outdir", default=".", help="results/ and logs/ go here")
    r.add_argument("--tmpdir", help="downloads and scratch (default OUTDIR/tmp)")
    r.add_argument("--index", help="bwa-mem3 index prefix (default: from the project file)")
    r.add_argument("--resolutions", help="comma list, or ladder:<base> for base*2^k up to 10 Mb (default: bin.resolutions)")
    r.add_argument("--filters", help="build only these of bin.filters (comma list)")
    r.add_argument("--cpus", type=int, help="threads to use (default: all)")
    r.add_argument("--mem", type=float, help="RAM budget in GB (default: 80%% of the machine)")
    r.add_argument("--downloads", type=int, default=3, help="parallel SRA downloads")
    r.add_argument("--sra-ahead", type=int, default=4, help="downloaded runs allowed to wait for the mapper")
    r.add_argument("--dump-procs", type=int, default=3, help="fastq-dump processes feeding the mapper")
    r.add_argument("--keep-sra", action="store_true", help="keep downloaded .sra files after mapping")
    r.add_argument("-n", "--dry-run", action="store_true", help="print the tasks and exit")
    pl = sub.add_parser("plan", help="print the project settings as JSON (used by main.nf)")
    for q in (r, pl)[1:]:
        q.add_argument("project"); q.add_argument("--index"); q.add_argument("--resolutions"); q.add_argument("--filters")
    s = sub.add_parser("sra-stream", help=".sra -> interleaved FASTQ on stdout")
    s.add_argument("sra"); s.add_argument("--procs", type=int, default=3)
    s.add_argument("--start", type=int); s.add_argument("--end", type=int)
    s.add_argument("--fastq-dump", default="fastq-dump"); s.add_argument("--vdb-dump", default="vdb-dump")
    t = sub.add_parser("sum-stats", help="sum pairtools-style stats files")
    t.add_argument("out"); t.add_argument("inputs", nargs="+")
    a = ap.parse_args()
    if a.cmd == "run": sys.exit(run_pipeline(a))
    if a.cmd == "plan":
        import json; print(json.dumps(plan(a), indent=1))
    if a.cmd == "sra-stream":
        try: sra_stream(a.sra, a.procs, a.start, a.end, a.fastq_dump, a.vdb_dump, sys.stdout.buffer)
        except BrokenPipeError: sys.exit(1)
    if a.cmd == "sum-stats": sum_stats(a.inputs, a.out)


if __name__ == "__main__":
    main()
