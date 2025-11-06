#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

Channel
    .fromList(['SRR000073', 'SRR000074', 'SRR000075'])
    .set { sra_ids_ch }

sudo apt update
sudo apt install sra-toolkit

process get_fastq {
    tag "$sra_id"
    publishDir "data"

    input:
    val sra_id

    output:
    path "${sra_id}.fastq.gz" into sra_ids_ch

    script:
    """
    fasterq-dump --threads 4 --progress ${sra_id}
    gzip *.fastq
    """
}

process trim_reads {
    tag { reads.simpleName }
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
    sra_ids_ch | get_fastq
    
}
