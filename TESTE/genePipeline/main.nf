#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.sra_ids = ['SRR10379721', 'SRR10379722', 'SRR10379723', 'SRR10379724', 'SRR10379725', 'SRR10379726']


Channel
	.fromList(params.sra_ids)
	.set { sra_ids_ch }


process download_and_gzip_fastq {
    tag { sra_id }
    container = "build('./Dockerfile.sratoolkit')"
    
    input:
    val sra_id from sra_ids_ch

    output:
    file "${sra_id}.fastq.gz" into fastq_ch

    script:
    """
    fasterq-dump --threads 1 --progress ${sra_id}
    gzip ${sra_id}.fastq
    """
}

process trim_reads {
    tag { reads.simpleName }
    container = "build('./Dockerfile.catadapt')"
    publishDir "results/trimming", mode: 'copy'

    input:
    file reads from fastq_ch

    output:
    file "${reads.simpleName}_trimmed.fastq.gz" into trimmed_ch

    script:
    """
    trim_galore -q 20 --phred33 --length 25 --output_dir . ${reads}
    mv ${reads.simpleName}_trimmed.fq.gz ${reads.simpleName}_trimmed.fastq.gz
    """
}


process get_reference {
    container = build('./Dockerfile.bowtie') 
    output:
    file 'reference.fasta'
    file 'reference.gff'
    script:
    """
    wget -O reference.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=fasta&retmode=text"
    wget -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=gff3&retmode=text"
    """
}


process build_index {
    container = build('./Dockerfile.bowtie')
    input:
    file fasta from get_reference.out.filter{ it.name.endsWith('.fasta') }
    output:
    file 'reference_index.*.ebwt'
    script:
    """
    bowtie-build reference.fasta reference_index
    """
}


workflow {
    download_and_gzip_fastq()
    trim_reads()
   # reference_files = get_reference()
   # index_files = buid_index(reference_files)
}

