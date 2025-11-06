#!/usr/bin/env nextflow

#params.sra_ids = ['SRR10379721', 'SRR10379722', 'SRR10379723', 'SRR10379724', 'SRR10379725', 'SRR10379726']
params.sra_ids = ['SRR000073']



process download_and_gzip_fastq {
    input:
    val sra_id from params.sra_ids

    output:
    file "${sra_id}.fastq.gz" into gzipped_files

    script:
    """
    fasterq-dump --threads 4 --progress ${sra_id}
    gzip ${sra_id}.fastq
    """
}

process trim_reads {
    input:
    val sra_id from params.sra_ids

    output:
    file "${sra_id}.fastq.gz" into trimmed_files

    container= build('./Dockerfile.catadapt')

    script:
    """
    trim_galore -q 20 --phred33 --length 25 ${sra_id}.fastq
    """
}

workflow {
    download_and_gzip_fastq()
    trim_reads()
}

