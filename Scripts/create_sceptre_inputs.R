#!/usr/bin/env Rscript
# Create SCEPTRE input files (matrix.mtx, features.tsv, barcodes.tsv)
# Snakemake inputs: guide_cdna, filtered_cells
# Snakemake outputs: matrix, features, barcodes

library(dplyr)
library(Matrix)
library(biomaRt)

# Read inputs from Snakemake
guide_cdna_file <- snakemake@input[['guide_cdna']]
filtered_cells_file <- snakemake@input[['filtered_cells']]

matrix_file <- snakemake@output[['matrix']]
features_file <- snakemake@output[['features']]
barcodes_file <- snakemake@output[['barcodes']]

gene_list <- snakemake@params[['gene_list']]

# Create output directory
dir.create(dirname(matrix_file), showWarnings = FALSE, recursive = TRUE)

# Read filtered data
input_file <- read.table(guide_cdna_file, header=TRUE, sep='\t')

# Get Ensembl IDs for genes using biomaRt
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id'),
  filters = 'hgnc_symbol',
  values = gene_list,
  mart = mart
)

gene_to_ensembl <- setNames(gene_map$ensembl_gene_id, gene_map$hgnc_symbol)

# Replace gene names with Ensembl IDs (except for guides)
input_file$gene <- sapply(input_file$gene, function(gene) {
  if (grepl("enh|_prom|DMC|safeTarg", gene)) {
    return(gene)  # keep original name if it's a guide
  }
  gene_base <- as.character(gene)
  new_name <- gene_to_ensembl[gene_base]
  if (is.na(new_name)) return(gene_base)
  return(new_name)
})

# Get unique genes and cells
genes_filtered <- unique(input_file$gene)
input_file$Cell_Barcode <- trimws(as.character(input_file$Cell_Barcode))
cells_filtered <- unique(input_file$Cell_Barcode)

# Create sparse matrix
gene_to_idx <- match(input_file$gene, genes_filtered)
cell_to_idx <- match(input_file$Cell_Barcode, cells_filtered)

matrix_filtered <- sparseMatrix(
  i = gene_to_idx,
  j = cell_to_idx,
  x = input_file$count,
  dims = c(length(genes_filtered), length(cells_filtered))
)

# Write matrix.mtx
writeMM(matrix_filtered, matrix_file)

# Build features.tsv with 3 columns
features_df <- data.frame(
  V1 = genes_filtered,
  V2 = genes_filtered,
  V3 = sapply(genes_filtered, function(gene) {
    if (grepl("enh|_prom|DMC|safeTarg", gene)) {
      return("CRISPR Guide Capture")
    }
    return("Gene Expression")
  }),
  stringsAsFactors = FALSE
)

# Replace ampersands with dashes (SCEPTRE requirement)
features_df$V1 <- gsub("&", "-", features_df$V1)
features_df$V2 <- gsub("&", "-", features_df$V2)

# Write output files
write.table(features_df, features_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
write.table(cells_filtered, barcodes_file, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
