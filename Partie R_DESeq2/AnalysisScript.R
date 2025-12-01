#!/usr/bin/env Rscript

###########################################################################################
# THIS PART MUST BE CHANGED IN ORDER TO PLACE INSIDE NEXT FLOW

# Diretório onde estão os arquivos baixados
# (ajusta conforme onde você salvou)
main_path <- getwd()   # ou: "/caminho/para/sua/pasta"

# Arquivos de entrada (baixados da máquina virtual)
count_table_path  <- file.path(main_path, "counts.txt")          # nome do seu counts
coldata_file_path <- file.path(main_path, "coldata.txt")         # nome do seu coldata
geneDB_path       <- file.path(main_path, "GeneSpecificInformation_NCTC8325.xlsx")

# Arquivos de saída (vão ser criados nesse diretório)
output_finalcountdata_path <- file.path(main_path, "deseq_input_countdata.csv")
output_vst_path            <- file.path(main_path, "vst_table.csv")
output_file_path           <- file.path(main_path, "deseq_results.csv")

############################################################################################


# Packages

library(DESeq2)
library(stringi)
library(stringr)
library(dplyr)
library(glue)
library(readxl)
library(ggplot2)
library(KEGGREST)
library(ggrepel)



## Functions

# This function arranges the count table for the analysis

formatation_table <- function(count_table, coldata) {
  columns    <- colnames(count_table)    # 
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

  # Connects each .bam column to your sample
  final_columns <- rep("", length(sample_columns))
  for (j in seq_along(sample_columns)) {
    for (i in seq_along(samples)) {
      if (isTRUE(str_detect(sample_columns[j], coldata$samples[i]))) {    #Detection of the samples
        final_columns[j] <- samples[i]
      }
    }
  }
  
  # correction of the columns names: sample_columns -> final_columns
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

# Function used to avoid repetition of gene names

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

coldata <- read.table(coldata_file_path, sep = " ", header = TRUE, check.names = FALSE)
coldata$condition <- as.character(coldata$condition)
coldata$samples   <- as.character(coldata$samples)

# Formatation of the count_table, and sync to coldata
formatted <- formatation_table(raw_count_table, coldata)

#Split formatted and update count_table (creation of a new table) and coldata
format_count_table <- formatted$countdata
coldata <- formatted$coldata


## From AureoWik: GeneDB (.xlsx)

# Reads file and replace "na" for "-"
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

#---------------------------------------------------------------------------------------------
# Merge counts table + gene names

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

## Salva matriz final de contagens usada pelo DESeq2
write.csv(final_count_table, file = output_finalcountdata_path, row.names = TRUE)

#--------------------------------------------------------------------------------------------------------
# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = final_count_table,
  colData   = coldata,
  design    = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)

#----------------------------------------------------------------------------------------------------
# MA plot + results

res <- results(dds, alpha = 0.05)
res$color <- ifelse(res$padj < 0.05, "red", "black")

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

pdf(file = file.path(main_path, "MA_plot_all_genes.pdf"),
    width = 8, height = 7, family = "Helvetica", useDingbats = FALSE)
print(p_ma)
dev.off()


## Export Results table
res <- as.data.frame(subset(res, select = -color))
res$conv_name <- gene_names
write.csv(res, file = output_file_path, row.names = TRUE)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#TEST

## ============================
## MA-plot só de genes de tradução (KEGG)
## ============================

message("Construindo MA-plot de genes relacionados à tradução (KEGG)...")

# Copia resultados para um data.frame
res_df <- res

# -----------------------------
# 1) Buscar genes de tradução no KEGG
#    sao03010, sao00970 = translation-related
#    sao03012 = translation factors
# -----------------------------

list_kegg <- c("sao03010", "sao00970")
all_gene  <- list(ID = list(), NAME = list(), AA_tRNA = list())

for (KEgg in list_kegg) {
  gene_keg <- keggGet(KEgg)[[1]]
  genes_vec <- gene_keg$GENE  # vetor do tipo: ID1, desc1, ID2, desc2, ...
  
  for (i in seq(1, length(genes_vec), by = 2)) {
    id   <- genes_vec[i]
    desc <- genes_vec[i + 1]
    
    all_gene$ID   <- append(all_gene$ID,   id)
    all_gene$NAME <- append(all_gene$NAME, desc)
    
    # Heurística: TRUE para genes que NÃO são pseudogene / tRNA / ribosomal
    is_aars <- str_detect(
      desc,
      stringr::str_c(
        "\\b(",
        stringr::str_c(c("pseudogene", "tRNA-", "ribosomal", "Ribosomal"),
                       collapse = "|"),
        ")\\b"
      ),
      negate = TRUE
    )
    
    all_gene$AA_tRNA <- append(all_gene$AA_tRNA, is_aars)
  }
}

# More genes 
all_gene$ID      <- append(all_gene$ID,   "SAOUHSC_01203")
all_gene$NAME    <- append(all_gene$NAME, "ribonuclease III")
all_gene$AA_tRNA <- append(all_gene$AA_tRNA, FALSE)

# -----------------------------
# 2) Translation factors (sao03012)
# -----------------------------

tmp_file <- tempfile(fileext = ".keg")
download.file("https://www.genome.jp/kegg-bin/download_htext?htext=sao03012", destfile = tmp_file, quiet = TRUE)

# Read lines from KEGG files 
sao03012_raw <- readLines(tmp_file)

# Fica só com linhas que tenham SAOUHSC_ (ids de genes)
sao03012_raw <- sao03012_raw[ stringr::str_detect(sao03012_raw, "SAOUHSC_") ]

# Extrai diretamente os IDs SAOUHSC_XXXXX
sao03012_id <- stringr::str_extract(sao03012_raw, "SAOUHSC_\\d+")

# Marca quais desses têm "tRNA" no texto (aminoacyl-tRNA-related)
sao03012_AA_tRNA <- stringr::str_detect(sao03012_raw, "tRNA")


# Adiciona translation factors na lista geral de genes
all_gene$ID      <- append(all_gene$ID,   unlist(sao03012_id))
all_gene$AA_tRNA <- append(all_gene$AA_tRNA, unlist(sao03012_AA_tRNA))
# (NAME não é usado mais pra frente, então ignoramos aqui)

# -----------------------------
# 3) Filtrar resultados DESeq2 só para esses genes de tradução
# -----------------------------

gene_ids_translation <- unlist(all_gene$ID)
aars_ids             <- unlist(all_gene$ID)[unlist(all_gene$AA_tRNA)]

# Mantém só genes cujo rowname (Geneid) está em gene_ids_translation
res_trans <- res_df[rownames(res_df) %in% gene_ids_translation, ]

if (nrow(res_trans) == 0) {
  warning("Nenhum gene de tradução (KEGG) encontrado nos resultados. MA-plot de tradução não será gerado.")
} else {
  # Marca quais são AA-tRNA synthetases
  res_trans$AA_tRNA <- rownames(res_trans) %in% aars_ids
}  
  # -----------------------------
  # 4) Preparar colunas para plot
  # -----------------------------
  
  # Log2 da média (como na figura do artigo)
  res_trans$log2BaseMean <- log2(res_trans$baseMean + 1)
  
  # Significância pelo padj
  res_trans$Significance <- ifelse(
    !is.na(res_trans$padj) & res_trans$padj <= 0.05,
    "Significant",
    "Non-Significant"
  )
  
  # Genes a rotular explicitamente
  genes_to_label <- c("frr", "infA", "tsf", "infC", "infB", "pth")
  res_trans$id_label <- ifelse(
    res_trans$conv_name %in% genes_to_label,
    res_trans$conv_name,
    ""
  )
  
  # -----------------------------
  # 5) MA-plot de tradução
  # -----------------------------
  
  p_trans <- ggplot(res_trans) +
    geom_point(aes(x= log2BaseMean,y= log2FoldChange,fill= Significance, color= AA_tRNA),
      size   = 3,
      shape  = 21,
      stroke = 0
    ) +
    scale_fill_manual(
      values = c(
        "Significant"     = "red",
        "Non-Significant" = "black"
      )
    ) +
    scale_color_manual(
      values = c(
        "TRUE"  = "black"  # AA-tRNA synthetases
      )
    ) +
    geom_text_repel(
      aes(
        x     = log2BaseMean,
        y     = log2FoldChange,
        label = id_label
      ),
      show.legend  = FALSE,
      size         = 6,
      box.padding  = 0.4,
      max.overlaps = 5000
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
  
  pdf(
    file   = file.path(main_path, "MA_plot_translation.pdf"),
    width  = 8,
    height = 7,
    family = "Helvetica",
    useDingbats = FALSE
  )
  print(p_trans)
  dev.off()
  
  message("MA_plot_translation.pdf salvo em: ", file.path(main_path, "MA_plot_translation.pdf"))




#
#
#
#
#


