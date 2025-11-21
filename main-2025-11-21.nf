#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
Channel
    .fromList(['SRR000073'])
    //.fromList(['SRR000073', 'SRR000074', 'SRR000075'])
    .set { sra_ids_ch }
process get_fastq {
    tag "$sra_id"
    publishDir "data/raw_fastq/", mode: 'copy'

    input:
    val sra_id

    output:
    path "${sra_id}.fastq.gz", emit: fastq_gz_ch

    script:
    """
    fasterq-dump --threads ${task.cpus} --progress ${sra_id}
    gzip *.fastq
    """
}
process trim_reads {
    publishDir "data/trimmed_fastq/", mode: 'copy'

    input:
    path fastq

    output:
    path "${fastq.baseName}_trimmed.fq.gz", emit: trimmed_ch
    path "cutadapt_version.txt"

    script:
    """
    cutadapt --version > cutadapt_version.txt
    trim_galore --cores ${task.cpus} -q 20 --phred33 --length 25 --basename ${fastq.baseName} ${fastq}
    """
}
process get_reference {
    publishDir "data/ref/", mode: 'copy'

    output:
    path 'reference.fasta', emit: reference_fasta_ch
    path 'reference.gff', emit: reference_gff_ch

    script:
    """
    wget -O reference.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=fasta&retmode=text"
    wget -O reference.gff "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=gff3&retmode=text"
    """
}
process build_index {
    publishDir "data/reference_index/", mode: 'copy'

    input:
    path fasta

    output:
    path 'reference_index.*.ebwt', emit: built_index_ch

    script:
    """
    bowtie-build ${fasta} reference_index
    """
}
process mapping {
    publishDir "data/map/", mode: 'copy'

    input:
    path built_index_basename
    path fastq_gz

    output:
    path "*.bam", emit: bam_ch

    script:
    """
    bowtie --threads ${task.cpus} -S ${built_index_basename} ${fastq_gz} | samtools sort --threads ${task.cpus} -o ${fastq_gz.baseName}.bam
    samtools index ${fastq_gz.baseName}.bam
    """
}
process feature_counts {
    publishDir "data/counts/", mode: 'copy'

    input:
    path reference_gff
    path bam_ch

    output:
    path 'counts.txt', emit: counts_ch

    script:
    """
    featureCounts -T ${task.cpus} -t gene -g ID -F GTF -a ${reference_gff} -o counts.txt ${bam_ch}
    """
}
workflow {
    get_fastq(sra_ids_ch)
    trim_reads(get_fastq.out.fastq_gz_ch)
    get_reference()
    build_index(get_reference.out.reference_fasta_ch)
    index_basename_ch = build_index.out.built_index_ch.map { list -> list[0].toString().split('\\.')[0] }
    mapping(index_basename_ch, get_fastq.out.fastq_gz_ch)
    mapping.out.bam_ch.collect() | set { all_bams_ch }
    feature_counts(get_reference.out.reference_gff_ch, all_bams_ch)
}
