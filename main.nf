#!/usr/bin/env nextflow
// rustiller on Nextflow: the same pipeline as `rustiller.py run`, for clusters.
//
//   nextflow run main.nf --project project.yml --outdir OUT [-profile slurm] [--resolutions ladder:1024]
//
// The project file is read by `rustiller.py plan`, so both executors interpret it identically.
// Paths in it (index, chromsizes, FASTQ) must be visible to the compute nodes.
nextflow.enable.dsl = 2

params.project = null
params.outdir = '.'
params.index = null
params.resolutions = null
params.filters = null
params.dump_procs = 3
params.publish_mode = 'copy'
params.python = 'python3'

def loadPlan() {
    if (!params.project) error "give --project <distiller project yaml>"
    def cmd = [params.python, "${projectDir}/rustiller.py", 'plan', file(params.project).toString()]
    if (params.index) cmd += ['--index', params.index.toString()]
    if (params.resolutions) cmd += ['--resolutions', params.resolutions.toString()]
    if (params.filters) cmd += ['--filters', params.filters.toString()]
    def p = cmd.execute()
    def out = new StringBuilder()
    def err = new StringBuilder()
    p.waitForProcessOutput(out, err)
    if (err) log.warn err.toString().trim()
    if (p.exitValue() != 0) error "rustiller.py plan failed"
    return new groovy.json.JsonSlurper().parseText(out.toString())
}

process FETCH_SRA {
    label 'download'
    tag "$acc"
    input: val acc
    output: tuple val(acc), path("${acc}.sra")
    script:
    """
    curl -sS -f --retry 8 --retry-delay 10 -o ${acc}.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/${acc}/${acc}
    """
}

process MAP_PARSE {
    label 'map'
    tag "${L.library}.${L.lane}"
    publishDir { "${params.outdir}/${P.dirs.mapped_parsed_sorted_chunks}" }, mode: params.publish_mode
    input: tuple val(L), path(reads), val(P)
    output: tuple val(L), path("${L.library}.${L.lane}.${P.assembly}.pairs.parquet")
    script:
    def bwa_t = Math.max(1, (task.cpus * 3).intdiv(4))
    def parse_t = Math.max(1, Math.min(4, task.cpus.intdiv(8)))
    def range = L.end ? "--start ${L.start} --end ${L.end}" : ''
    def feed = L.sra ? "rustiller.py sra-stream ${reads} --procs ${params.dump_procs} --fastq-dump '${P.tools['fastq-dump']}' --vdb-dump '${P.tools['vdb-dump']}' ${range} |" : ''
    def bwa_in = L.sra ? '-p -' : reads.join(' ')
    def out = "${L.library}.${L.lane}.${P.assembly}.pairs.parquet"
    """
    set -o pipefail
    ${feed} ${P.tools['bwa-mem3']} mem -SP5M ${P.map_opts} -t ${bwa_t} ${P.index} ${bwa_in} \\
      | ${P.tools.rooler} pairs parse - part.parquet --chroms-path ${P.chromsizes} --assembly ${P.assembly} ${P.parse_opts} --threads ${parse_t} --mem 8
    mv part.parquet ${out}
    """
}

process DEDUP {
    label 'helper'
    tag "$lib"
    publishDir { "${params.outdir}/${P.dirs.pairs_library}" }, mode: params.publish_mode, pattern: '*.parquet'
    publishDir { "${params.outdir}/${P.dirs.stats_library}" }, mode: params.publish_mode, pattern: '*.stats'
    input: tuple val(lib), path(lanes), val(P)
    output:
    tuple val(lib), path("${lib}.${P.assembly}.nodups.pairs.parquet"), emit: pairs
    tuple val(lib), path("${lib}.${P.assembly}.dedup.stats"), emit: stats
    script:
    """
    ${P.tools.rooler} pairs merge ${lib}.${P.assembly}.nodups.pairs.parquet ${lanes} --dedup ${P.max_mismatch} \\
      --stats ${lib}.${P.assembly}.dedup.stats --threads ${task.cpus} --mem ${Math.max(2, task.memory.toGiga() - 4)}
    """
}

process CLOAD {
    label 'helper'
    tag "${lib}.${fname}"
    publishDir { "${params.outdir}/${P.dirs.coolers_library}" }, mode: params.publish_mode
    input: tuple val(lib), path(pairs), val(fname), val(minq), val(P)
    output: tuple val(lib), val(fname), path("${lib}.${P.assembly}.${fname}.${P.resolutions[0]}.cool")
    script:
    """
    set -o pipefail
    ${P.tools.rooler} pairs view ${pairs} --min-mapq ${minq} --chroms ${P.chromsizes} \\
      | ${P.tools.rooler} cload - ${P.resolutions[0]} ${lib}.${P.assembly}.${fname}.${P.resolutions[0]}.cool \\
          --assembly ${P.assembly} --threads ${task.cpus} --mem ${Math.max(2, task.memory.toGiga() - 4)}
    """
}

process MERGE_COOLERS {
    label 'helper'
    tag "${group}.${fname}"
    publishDir { "${params.outdir}/${P.dirs.coolers_library_group}" }, mode: params.publish_mode
    input: tuple val(group), val(fname), path(cools), val(P)
    output: tuple val(group), val(fname), path("${group}.${P.assembly}.${fname}.${P.resolutions[0]}.cool")
    script:
    """
    ${P.tools.rooler} merge ${group}.${P.assembly}.${fname}.${P.resolutions[0]}.cool ${cools} \\
      --assembly ${P.assembly} --threads ${task.cpus} --mem ${Math.max(2, task.memory.toGiga() - 4)}
    """
}

process ZOOMIFY {
    label 'helper'
    tag "${cool.baseName}"
    publishDir { "${params.outdir}/${P.dirs[dirkey]}" }, mode: params.publish_mode
    input: tuple val(dirkey), path(cool), val(P)
    output: path("${cool.baseName}.mcool")
    script:
    """
    ${P.tools.rooler} zoomify ${cool} ${cool.baseName}.mcool --resolutions ${P.resolutions.join(',')} ${P.balance ? '--balance' : ''} \\
      --assembly ${P.assembly} --threads ${task.cpus} --mem ${Math.max(2, task.memory.toGiga() - 4)}
    """
}

process GROUP_STATS {
    tag "$group"
    publishDir { "${params.outdir}/${P.dirs.stats_library_group}" }, mode: params.publish_mode
    input: tuple val(group), path(stats), val(P)
    output: path("${group}.${P.assembly}.stats")
    script:
    """
    rustiller.py sum-stats ${group}.${P.assembly}.stats ${stats}
    """
}

workflow {
    P = loadPlan()
    nLanes = P.lanes.countBy { it.library }
    groupsOf = [:].withDefault { [] }                      // library -> groups it belongs to
    P.groups.each { g, libs -> libs.each { groupsOf[it] << g } }

    lanes = Channel.fromList(P.lanes)
    sra = FETCH_SRA(lanes.filter { it.sra }.map { it.sra }.unique())
    reads = lanes.filter { it.sra }.map { [it.sra, it] }.combine(sra, by: 0).map { acc, L, f -> [L, f] }
        .mix(lanes.filter { !it.sra }.map { [it, it.fastqs.collect { f -> file(f) }] })
    mapped = MAP_PARSE(reads.map { L, f -> [L, f, P] })

    // groupKey carries the expected size, so a library moves on when its own lanes are done
    // instead of waiting for every lane of the project
    DEDUP(mapped.map { L, pq -> [groupKey(L.library, nLanes[L.library]), pq] }.groupTuple()
                .map { lib, pqs -> [lib.toString(), pqs, P] })

    filters = Channel.fromList(P.filters.collect { name, q -> [name, q] })
    cools = CLOAD(DEDUP.out.pairs.combine(filters).map { lib, pq, fname, q -> [lib, pq, fname, q, P] })
    gcools = MERGE_COOLERS(
        cools.flatMap { lib, fname, cool -> groupsOf[lib].collect { g -> [groupKey("${g}\t${fname}".toString(), P.groups[g].size()), cool] } }
             .groupTuple()
             .map { key, cs -> def (g, fname) = key.toString().split('\t'); [g, fname, cs, P] })
    ZOOMIFY(cools.map { lib, fname, cool -> ['coolers_library', cool, P] }
                 .mix(gcools.map { g, fname, cool -> ['coolers_library_group', cool, P] }))

    GROUP_STATS(DEDUP.out.stats.flatMap { lib, st -> groupsOf[lib].collect { g -> [groupKey(g, P.groups[g].size()), st] } }
                     .groupTuple().map { g, sts -> [g.toString(), sts, P] })
}
