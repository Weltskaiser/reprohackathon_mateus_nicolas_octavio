#!/usr/bin/env Rscript

##################################################
## Debug on personal computer, without Nextflow ##
####################################################################################
# Working directory, containing useful data
# If necessary, change according to where they are
#main_path <- getwd()   # ou: "/caminho/para/sua/pasta"

# Input files
#count_table_path  <- file.path(main_path, "counts.txt")
#geneDB_path       <- file.path(main_path, "GeneSpecificInformation_NCTC8325.xlsx")
# KEGG Orthology
#kegg_path         <- file.path(main_path, "kegg_2_sao.json")
# Output files (will be created in the same directory)
#output_finalcountdata_path <- file.path(main_path, "deseq_input_countdata.csv")
#output_deseq_results_path           <- file.path(main_path, "deseq_results.csv")
####################################################################################

args <- commandArgs(trailingOnly = TRUE)
# Input paths
count_table_path              <- args[1]
geneDB_path                   <- args[2]
kegg_path                     <- args[3]
# Output paths
output_finalcountdata_path    <- args[4]
output_deseq_results_path     <- args[5]
# Working directory
main_path                     <- args[6]


# Packages loading
library(DESeq2)
library(stringi)
library(stringr)
library(dplyr)
library(glue)
library(readxl)
library(ggplot2)
library(ggrepel)
library(jsonlite)


# Arrange the count table for the analysis
formatation_table <- function(count_table, coldata) {
  columns    <- colnames(count_table)
  conditions <- unique(coldata$condition)    # Control or  treatment
  samples    <- rep("", length(coldata$condition))

  # Create the names for each sample: control_1..., treatment_1...
  # Concatenates condition + _ + iteration
  iter <- 1
  for (cond in seq_along(conditions)) {
    for (i in seq_along(coldata$condition)) {
      sample_condition <- coldata$condition[i]
      if (sample_condition == conditions[cond]) {
        samples[i] <- paste0(sample_condition, "_", iter)
        iter <- iter + 1
      }
    }
  }

  # Detection of the samples columns from counts.txt (now columns)
  sample_columns <- columns[str_detect(columns, ".bam")]

  # Connects each .bam column to the sample
  final_columns <- rep("", length(sample_columns))
  for (j in seq_along(sample_columns)) {
    for (i in seq_along(samples)) {
      if (isTRUE(str_detect(sample_columns[j], coldata$samples[i]))) {    #Detection of the samples
        final_columns[j] <- samples[i]
      }
    }
  }

  # Correction of the columns names: sample_columns -> final_columns
  table <- count_table[, c("Geneid", sample_columns)]
  names(table) <- c("Geneid", final_columns)
  table$Geneid <- str_replace(table$Geneid, "gene-", "")   # Erases "gene" from the sample names

  cols <- c("Geneid", samples)
  table <- table[, c(cols, setdiff(names(table), cols))]

  # Update samples, output countdata + coldata
  coldata$samples <- samples
  output <- list("countdata" = table, "coldata" = coldata)
  return(output)
}

# Avoid repetition of gene names
avoid_repetition_name <- function(data, column) {
  column <- deparse(substitute(column))
  data$occurrence <- ave(seq_len(nrow(data)), data[[column]], FUN = seq_along)

  data[[column]] <- ifelse(
    data$occurrence == 1,
    data[[column]],
    paste0(data[[column]], "_", data$occurrence - 1)
  )

  data$occurrence <- NULL
  return(data)
}


# Load tables
raw_count_table <- read.table(count_table_path, sep = "\t", header = TRUE, check.names = FALSE)

sample_id_l = c("SRR10379721", "SRR10379722", "SRR10379723", "SRR10379724", "SRR10379725", "SRR10379726")
condition_l = c("treatment", "treatment", "treatment", "control", "control", "control")
coldata <- data.frame(
  samples   = sample_id_l,
  condition = condition_l,
  stringsAsFactors = FALSE
)

# Formatation of the count_table, and sync to coldata
formatted <- formatation_table(raw_count_table, coldata)

# Split formatted and update count_table (creation of a new table) and coldata
format_count_table <- formatted$countdata
coldata <- formatted$coldata


## From AureoWik: GeneDB (.xlsx)

# Read file and replace "na" for "-"
ID_to_Names <- read_excel(geneDB_path)
ID_to_Names$`pan gene symbol`[is.na(ID_to_Names$`pan gene symbol`)] <- "-"

# If there is a name, name; if not, locus tag
gene_names <- rep("", length(ID_to_Names$`pan gene symbol`))
for (row in seq_along(ID_to_Names$`pan gene symbol`)) {
  if (ID_to_Names$`pan gene symbol`[row] == "-") {
    gene_names[row] <- ID_to_Names$`locus tag`[row]
  } else {
    gene_names[row] <- ID_to_Names$`pan gene symbol`[row]
  }
}

#Comment
ID_to_Names$Name <- gene_names
ID_to_Names <- ID_to_Names %>% dplyr::select(`locus tag`, Name)
colnames(ID_to_Names) <- c("Geneid", "Name")
ID_to_Names <- avoid_repetition_name(ID_to_Names, Name)


## Merge counts table + gene names

# Join through "Geneid"
merged <- inner_join(ID_to_Names, format_count_table, by = "Geneid")
final_count_table <- merged

# Remove gens where counts = 0
final_count_table <- final_count_table[
  rowSums(final_count_table[, !colnames(final_count_table) %in% c("Name", "Geneid")]) > 0,
]

row_names  <- final_count_table$Geneid

# This is used later on for Bio Annotation
gene_names <- final_count_table$Name

final_count_table <- final_count_table[, !names(final_count_table) %in% c("Name", "Geneid")]
rownames(final_count_table) <- row_names

## Export the final counts matrix used by DESeq2
write.csv(final_count_table, file = output_finalcountdata_path, row.names = TRUE)



## DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = final_count_table,
  colData   = coldata,
  design    = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)

res <- results(dds, alpha = 0.05)
res$color <- ifelse(res$padj < 0.05, "red", "black")


## 1st MA plot with all genes + results

p_ma <- ggplot(as.data.frame(res), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_manual(values = c("red" = "red", "black" = "black")) +
  theme_grey() +
  labs(
    title = "MA Plot of Differentially Expressed Genes",
    x = "Mean of normalized counts",
    y = "Log2 Fold Change"
  ) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  coord_cartesian(xlim = c(0.1, 4 * 10^5), ylim = c(-4, 4)) +
  scale_x_log10()

# Export as PDF
pdf(file = file.path(main_path, "MA_plot_all_genes.pdf"),
    width = 8, height = 7, family = "Helvetica", useDingbats = FALSE)
print(p_ma)
dev.off()

# Export Results table
res <- as.data.frame(subset(res, select = -color))
res$conv_name <- gene_names
write.csv(res, file = output_deseq_results_path, row.names = TRUE)




## MA plot with translation genes

# Make a copy of the results from DESeq2
res_df <- res

# Search translation genes in KEGG Orthology
# sao03010, sao00970 = translation-related
# sao03012 = translation factors

# We made an archive of the useful data from KEGG, in case it would change
# Refer to get_kegg_2_sao.R to see how it has been built
kegg_2_sao <- fromJSON(kegg_path)

list_kegg <- c("sao03010", "sao00970")
all_gene  <- list(ID = character(), NAME = character(), AA_tRNA = logical())

for (kegg in list_kegg) {
  gene_keg <- kegg_2_sao[[kegg]]
  genes_vec <- gene_keg$GENE[[1]]  # Vector of type: ID1, desc1, ID2, desc2, ...

  for (i in seq(1, length(genes_vec), by = 2)) {
    id   <- genes_vec[i]
    desc <- genes_vec[i + 1]

    all_gene$ID   <- c(all_gene$ID,   id)
    all_gene$NAME <- c(all_gene$NAME, desc)

    # Heuristic: TRUE for genes that are not pseudogenes / tRNA / ribosomal
    is_aars <- !str_detect(
      desc,
      stringr::str_c(
        "\\b(",
        stringr::str_c(c("-tRNA"),
                       collapse = "|"),
        ")\\b"
      )
    )

    all_gene$AA_tRNA <- c(all_gene$AA_tRNA, is_aars)
  }
}

# More genes
all_gene$ID      <- append(all_gene$ID,   "SAOUHSC_01203")
all_gene$NAME    <- append(all_gene$NAME, "ribonuclease III")
all_gene$AA_tRNA <- append(all_gene$AA_tRNA, FALSE)

# Translation factors (sao03012)

tmp_file <- tempfile(fileext = ".keg")
download.file("https://www.genome.jp/kegg-bin/download_htext?htext=sao03012", destfile = tmp_file, quiet = TRUE)

# Read data from KEGG
sao03012_raw <- readLines(tmp_file)
# Keep only lines that have SAOUHSC_ (gene IDs)
sao03012_raw <- sao03012_raw[ stringr::str_detect(sao03012_raw, "SAOUHSC_") ]
# Extract SAOUHSC_XXXXX IDs directly
sao03012_id <- stringr::str_extract(sao03012_raw, "SAOUHSC_\\d+")
# Mark which ones have "tRNA" in the text (aminoacyl-tRNA-related)
sao03012_AA_tRNA <- stringr::str_detect(sao03012_raw, "tRNA")

# Add translation factors to the general gene list
all_gene$ID      <- append(all_gene$ID,   unlist(sao03012_id))
all_gene$AA_tRNA <- append(all_gene$AA_tRNA, unlist(sao03012_AA_tRNA))
# (NAME is not used later, so we ignore it here)

# Filter DESeq2 results only for these translation genes

gene_ids_translation <- unlist(all_gene$ID)
aars_ids             <- unlist(all_gene$ID)[unlist(all_gene$AA_tRNA)]

# Keep only genes whose rowname (Geneid) is in gene_ids_translation
res_trans <- res_df[rownames(res_df) %in% gene_ids_translation, ]

if (nrow(res_trans) == 0) {
  warning("Nenhum gene de tradução (KEGG) encontrado nos resultados. MA-plot de tradução não será gerado.")
} else {
  # Mark which ones are AA-tRNA synthetases
  res_trans$AA_tRNA <- rownames(res_trans) %in% aars_ids
}

# Prepare columns for plot

# Log2 of the mean (as in the figure in the article)
res_trans$log2BaseMean <- log2(res_trans$baseMean + 1)

# Significance by padj
res_trans$Significance <- ifelse(
  !is.na(res_trans$padj) & res_trans$padj <= 0.05,
  "Significant",
  "Non-Significant"
)

# Genes to explicitly label
genes_to_label <- c("frr", "infA", "tsf", "infC", "infB", "pth")
res_trans$id_label <- ifelse(
  res_trans$conv_name %in% genes_to_label,
  res_trans$conv_name,
  ""
)


## MA plot with translation genes

p_trans <- ggplot(res_trans) +
  geom_point(aes(x= log2BaseMean,y= log2FoldChange,fill= Significance, color= AA_tRNA),
    size   = 3,
    shape  = 21,
    stroke = 1
  ) +
  scale_fill_manual(
    values = c(
      "Significant"     = "red",
      "Non-Significant" = "gray"
    )
  ) +
  scale_color_manual(
    values = c(
      "TRUE"  = "black",  # AA-tRNA synthetases
      "FALSE" = "NA"      # Non-AA-tRNA synthetases
    ),
    labels = c("TRUE" = "AA-tRNA synthethases"),
    breaks = "TRUE"
  ) +
  geom_text_repel(
    aes(
      x     = log2BaseMean,
      y     = log2FoldChange,
      label = id_label
    ),
    show.legend  = FALSE,
    size         = 6,
    box.padding  = 0.4
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(breaks = seq(-5, 5, 1)) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  labs(
    x = "Log2 base Mean",
    y = "Log2 Fold Change",
    title = "MA-plot of translation-related genes (KEGG)"
  ) +
  theme_minimal()

# Export as PDF
pdf(
  file   = file.path(main_path, "MA_plot_translation.pdf"),
  width  = 8,
  height = 7,
  family = "Helvetica",
  useDingbats = FALSE
)
print(p_trans)
dev.off()
