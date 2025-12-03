#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

Channel
    .fromList(['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726'])
    //.fromList(['SRR10379721','SRR10379722','SRR10379723','SRR10379724','SRR10379725','SRR10379726','SRR000073', 'SRR000074', 'SRR000075'])
    .set { sra_ids_ch }

Channel
    .fromPath('GeneSpecificInformation_NCTC8325.xlsx')
    .set { geneDB_ch }

Channel
    .fromPath('differential_analysis/differential_analysis.R')
    .set { analysis_script_ch }


def CONTROL_IDS = ['SRR10379724','SRR10379725','SRR10379726']


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
    gzip ${sra_id}.fastq
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
    tuple val("reference_index"), path("reference_index.*.ebwt"), emit: index_ch

    script:
    """
    bowtie-build ${fasta} reference_index
    """
}
process mapping {
    publishDir "data/map/", mode: 'copy'

    input:
    tuple val(index_name), path(index_files)
    path fastq_gz

    output:
    path "${fastq_gz.simpleName}.sorted.bam", emit: bam_ch

    script:
    """
    gunzip -c ${fastq_gz} > ${fastq_gz.simpleName}.fastq

    bowtie -p ${task.cpus} -S ${index_name} ${fastq_gz.simpleName}.fastq \
      | samtools view -bS - \
      | samtools sort -o ${fastq_gz.simpleName}.sorted.bam -
    samtools index ${fastq_gz.simpleName}.sorted.bam
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

process create_control_table {

    publishDir "data/coldata/", mode: 'copy'

    input:
    path bam_files      // lista de BAMs (vem de all_bams_ch)
    val  control_ids    // vetor com SRRs de controle (vem de CONTROL_IDS)

    output:
    path "coldata.txt", emit: coldata_ch

    script:
    """
    Rscript - << 'EOF'
    args_bam  <- "${bam_files.join(',')}"
    args_ctrl <- "${control_ids.join(',')}"

    mapped_files_path <- unlist(strsplit(args_bam, ","))
    control_ids        <- unlist(strsplit(args_ctrl, ","))

    extract_sra_id <- function(file_name) {
      match <- regmatches(file_name, gregexpr("SRR[0-9]+", file_name))
      if (length(match[[1]]) > 0) {
        return(match[[1]][1])
      } else {
        return(NA)
      }
    }

    sample_ids <- sapply(mapped_files_path, extract_sra_id)

    coldata <- data.frame(
      samples   = sample_ids,
      condition = ifelse(sample_ids %in% control_ids, "control", "treatment"),
      stringsAsFactors = FALSE
    )

    write.table(
      coldata,
      file      = "coldata.txt",
      sep       = " ",
      row.names = FALSE,
      quote     = FALSE
    )
    EOF
    """
}


process stat_analysis {

    publishDir "results/deseq2", mode: 'copy', overwrite: true

    input:
    path count_table
    path coldata_file
    path geneDB_file
    path analysis_script

    output:
    path "deseq_input_countdata.csv"
    path "vst_table.csv"
    path "deseq_results.csv"
    path "*.pdf"

    script:
    """
    Rscript ${analysis_script} \
      "${count_table}" \
      "${coldata_file}" \
      "deseq_input_countdata.csv" \
      "vst_table.csv" \
      "deseq_results.csv" \
      "${geneDB_file}" \
      \$PWD

    """
}


workflow {
    get_fastq(sra_ids_ch)
    trim_reads(get_fastq.out.fastq_gz_ch)
    get_reference()
    build_index(get_reference.out.reference_fasta_ch)
    index_ch = build_index.out.index_ch
    mapping(build_index.out.index_ch, trim_reads.out.trimmed_ch)
    mapping.out.bam_ch.collect() | set { all_bams_ch }
    feature_counts(get_reference.out.reference_gff_ch, all_bams_ch)
    create_control_table(all_bams_ch, CONTROL_IDS)
    stat_analysis(feature_counts.out.counts_ch, create_control_table.out.coldata_ch, geneDB_ch, analysis_script_ch)
}
