#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

Channel
    .fromList(['SRR000073', 'SRR000074', 'SRR000075'])
    .set { sra_ids_ch }

process get_fastq {
    tag "$sra_id"
    publishDir "data/raw_fastq/", mode: 'copy'

    input:
    val sra_id

    output:
    path "${sra_id}.fastq.gz", emit: fastq_files_ch

    script:
    """
    fasterq-dump --threads 4 --progress ${sra_id}
    gzip *.fastq
    """
}

process trim_reads {
    publishDir "data/trimmed_fastq/", mode: 'copy'

    input:
    path fastq_file

    output:
    path "${fastq_file.baseName}_trimmed.fq.gz", emit: trimmed_files
    path "cutadapt_version.txt", emit: version_txt

    container 'cutadapt:1.11'
    containerEngine = 'docker'

    script:
    """
    cutadapt --version > cutadapt_version.txt
    trim_galore -q 20 --phred33 --length 25 --basename ${fastq_file.baseName} ${fastq_file}
    """
}

workflow {
    get_fastq(sra_ids_ch)
    trim_reads(get_fastq.out.fastq_files_ch)
}
