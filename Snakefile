#!/usr/bin/env python3

import pandas

###########
# GLOBALS #
###########

run_info_file = 'data/SraRunInfo.csv'
csd_locus_fasta = 'data/csd.fasta'
bbduk_ref = '/sequencing_artifacts.fa.gz'
bbduk_adaptors = '/adapters.fa'
pipeline_threads = 20

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
gatk_container = 'shub://TomHarrop/singularity-containers:gatk_4.1.0.0'

########
# MAIN #
########

# get a list of all names
run_info = pandas.read_csv(run_info_file)
all_samples = sorted(set(run_info['SampleName']))

# TEMPORARY # subset samples to get the pipeline to run quickly
all_samples = [x for x in all_samples if x.startswith('REUd1')]

#########
# RULES #
#########

rule target:
    input:
        expand('output/030_markdup/{sample}.unsorted.duplicates_marked.bam',
               sample=all_samples)

# 030 gatk mark duplicates
rule markdup:
    input:
        'output/020_bwa/{sample}.sam'
    output:
        bam = 'output/030_markdup/{sample}.unsorted.duplicates_marked.bam',
        metrics = 'output/030_markdup/{sample}.duplicate_metrics'
    threads:
        pipeline_threads
    log:
        'output/000_logs/030_markdup/{sample}.log'
    benchmark:
        'output/000_benchmarks/030_markdup/{sample}.tsv'
    singularity:
        gatk_container
    shell:
        'gatk '
        'MarkDuplicates '
        '--INPUT {input} '
        '--OUTPUT {output.bam} '
        '--METRICS_FILE {output.metrics} '
        '--VALIDATION_STRINGENCY SILENT '
        '--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 '
        '--ASSUME_SORT_ORDER "queryname" '
        '--CREATE_MD5_FILE true '
        '&> {log}'


# 020 map each individual to reference contig
rule bwa:
    input:
        fq = 'data/reads/{sample}.fq.gz',
        index = expand('output/020_bwa/csd.fasta.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        'output/020_bwa/{sample}.sam'
    params:
        prefix = 'output/020_bwa/csd.fasta',
        rg = '\'@RG\\tID:{sample}\\tSM:{sample}\''
    threads:
        pipeline_threads
    log:
        'output/000_logs/020_bwa/{sample}.log'
    benchmark:
        'output/000_benchmarks/020_bwa/{sample}.tsv'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-K 100000000 '
        '-t {threads} '
        '-p '
        '-M '
        '-R {params.rg} '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

rule index:
    input:
        csd_locus_fasta
    output:
        expand('output/020_bwa/csd.fasta.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/020_bwa/csd.fasta'
    threads:
        1
    log:
        'output/000_logs/020_bwa/index.log'
    benchmark:
        'output/000_benchmarks/020_bwa/index.tsv'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'
