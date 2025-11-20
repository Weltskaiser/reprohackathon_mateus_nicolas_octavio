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

    container 'cutadapt'    
    containerEngine = 'docker'

    script:
    """
    cutadapt --version > cutadapt_version.txt
    trim_galore -q 20 --phred33 --length 25 --basename ${fastq_file.baseName} ${fastq_file}
    """
}


process get_reference {

    // se o wget está no host, você pode tirar o container;
    // se estiver dentro da imagem do bowtie, pode usar:
    // container 'bowtie2'

    output:
    path 'reference.fasta', emit: fasta
    path 'reference.gff',   emit: gff

    script:
    """
    wget -O reference.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=fasta&retmode=text"
    wget -O reference.gff   "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=CP000253.1&report=gff3&retmode=text"
    """
}

process build_index {

    container 'bowtie'   // nome da imagem que você construiu
    cpus 1

    input:
    path fasta
    path gff   // não é usado no comando, mas garante que o arquivo existe

    output:
    // todos os arquivos que o bowtie-build gera
    path 'reference_index.1.ebwt',    emit: idx1
    path 'reference_index.2.ebwt',    emit: idx2
    path 'reference_index.3.ebwt',    emit: idx3
    path 'reference_index.4.ebwt',    emit: idx4
    path 'reference_index.rev.1.ebwt', emit: rev1
    path 'reference_index.rev.2.ebwt', emit: rev2

    script:
    """
    bowtie-build ${fasta} reference_index
    """
}

//
// WORKFLOW PRINCIPAL
//

workflow {

    get_fastq(sra_ids_ch)
    trim_reads(get_fastq.out.fastq_files_ch)

    // baixa fasta + gff
    ref = get_reference()
    fasta_ch = ref.fasta
    gff_ch   = ref.gff

    // constrói o índice
    idx = build_index(fasta_ch, gff_ch)

    // prints só para ver que rodou
    fasta_ch.view { "FASTA: ${it}" }
    gff_ch.view   { "GFF:   ${it}" }
}

