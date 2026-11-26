# rustiller

A small, eager replacement for [distiller](https://github.com/open2c/distiller-nf): it reads a
distiller project file and produces deduplicated pairs, stats and balanced multi-resolution
coolers, using **bwa-mem3** for mapping and **rooler** for everything after it.

```
.sra ──fastq-dump (streamed)──▶ bwa-mem3 ──▶ rooler pairs parse ──▶ lane pairs (Parquet)
lanes ──rooler pairs merge --dedup──▶ library pairs + stats
library pairs ──rooler pairs view | rooler cload──▶ .cool ──rooler zoomify --balance──▶ .mcool
library coolers ──rooler merge──▶ group .cool ──▶ .mcool
```

## Two executors, one pipeline

| | `rustiller.py run` | `nextflow run main.nf` |
|---|---|---|
| For | one machine | a cluster (or anything Nextflow can submit to) |
| Scheduling | one mapper with most threads; helpers never delay it | one job per task, sized in `nextflow.config` |
| Files | written once, in place; downloads on `--tmpdir`, deleted after mapping | staged in `work/`, then published (copied) to the results |
| Needs | Python + PyYAML | also Nextflow (Java 17+) |

Both read the project file through the same code (`rustiller.py plan`), run the same commands
and give byte-identical pairs, stats and coolers. Both are eager.

## Use

```
rustiller.py run project.yml -o OUTDIR [--tmpdir FAST_DISK] [--resolutions ladder:1024] [--filters mapq_30]
rustiller.py run project.yml -o OUTDIR -n        # print the tasks and their commands

nextflow run main.nf --project project.yml --outdir OUTDIR [-profile slurm] [--resolutions ladder:1024] [-resume]
```

With Nextflow, the paths in the project file (index, chromsizes, FASTQ) must be visible to the
compute nodes, and `--publish_mode link` avoids the copy when `work/` and the results share a
filesystem. Task sizes (32 threads for the mapper, 8 for the rest) are in `nextflow.config`.

Re-running the same command resumes: a task whose output exists is skipped.

Needs Python with PyYAML, and on `PATH`: `rooler`, `bwa-mem3`, `fastq-dump`, `vdb-dump`
(sra-tools), `curl`. A tool can be pointed elsewhere in the project file or with an environment
variable (`RUSTILLER_ROOLER`, `RUSTILLER_BWA_MEM3`, ...):

```yaml
input:
    genome:
        bwa_mem3_index_path: /path/to/mm10.fa      # prefix of a `bwa-mem3 index` (not a bwa index)
rustiller:
    tools:
        rooler: /path/to/rooler
```

## What differs from distiller

- **Eager.** Each task starts when its own inputs are ready. A library is deduplicated and binned
  while later lanes map; distiller-nf waits for all mapping to finish. (In `main.nf` this is
  `groupKey(library, n_lanes)` in front of `groupTuple()`.)
- **No FASTQ on disk.** An SRA run is downloaded (AWS open-data mirror), then `fastq-dump`
  processes stream it straight into the mapper. Three of them deliver about 1.1 M read pairs/s;
  the mapper takes about 0.55 M/s on 32 threads.
- **One mapper at a time, with most of the threads.** On a 32-thread machine, one bwa-mem3
  process with 24 threads maps as fast as 2×16 or 3×10 and faster than 4×8, so lanes are mapped
  one after another and the remaining threads go to parsing and to one helper task (dedup,
  cload, zoomify). Helpers never hold threads the mapper is waiting for; once mapping is over
  they run in parallel.
- **Pairs are Parquet** (`*.pairs.parquet`, sorted, about 6.6 bytes per pair), not `.pairs.gz`.
  `rooler pairs view` turns one back into 4DN text.

## Project file support

| Key | Status |
|---|---|
| `input.raw_reads_paths` | `sra:<run>`, `sra:<run>?start=N&end=M`, or two local FASTQ paths per lane |
| `input.library_groups` | supported |
| `map.mapping_options` | passed to `bwa-mem3 mem -SP5M` |
| `map.chunksize`, `do_fastqc` | ignored |
| `map.trim_options` | not supported |
| `parse.parsing_options` | `--walks-policy all --add-columns mapq` is what rooler implements; `--min-mapq`, `--max-inter-align-gap`, `--max-molecule-size` are honoured; a warning is printed for anything else |
| `dedup.max_mismatch_bp` | supported |
| `bin.resolutions`, `bin.balance` | supported; `--resolutions` overrides |
| `bin.filters` | `''` and `'(mapq1>=N) and (mapq2>=N)'` |
| `output.dirs` | honoured for pairs, coolers and group stats |

## Notes

- bwa-mem3 breaks ties between equally good hits using the read's position in the input, so the
  order in which reads arrive changes about one pair per million. The stream is deterministic for
  a given `--dump-procs`, so a rerun gives identical output.
- Stats follow `pairtools stats` but have no `dist_freq` section.
