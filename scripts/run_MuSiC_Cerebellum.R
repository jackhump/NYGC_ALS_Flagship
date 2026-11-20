# Estimate cerebellum cell-type proportions with MuSiC
#
# Run from a working directory containing:
#   ene_matrix.RData
#   metadata.csv
#   GSE246183_230502_human_ctype.h5ad
#
# Required Bioconductor packages: MuSiC, biomaRt, zellkonverter, and
# SummarizedExperiment. Install compatible versions and dependencies, e.g.:
#   install.packages("BiocManager")
#   BiocManager::install(c("biomaRt", "zellkonverter", "SummarizedExperiment"))
# Install MuSiC using its project installation instructions.
#
# Input requirements:
# - The RData file contains genes_counts, a data frame with versioned Ensembl
#   gene IDs as row names and bulk RNA sample IDs as column names.
# - The metadata CSV has sample IDs in its first column and a Sample.Source
#   column. All count samples are assumed to have matching metadata.
# - The H5AD reference contains genotype, batch, and sub4 cell metadata and
#   an assay named X suitable for use as MuSiC counts.
#
# The script selects cerebellum samples, maps bulk Ensembl IDs to HGNC symbols
# using Ensembl, and estimates cell-type proportions using control reference
# cells. The Ensembl query requires an internet connection; mappings can change
# as Ensembl is updated.
#
# A converted reference RDS and the weighted proportion matrix CSV are written
# to the working directory. Existing files with the same names are overwritten.

required_packages <- c(
  "MuSiC", "biomaRt", "zellkonverter", "SummarizedExperiment"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

# Load counts separately so unrelated objects in the RData file do not enter
# the main workspace.
count_data <- new.env()
load("ene_matrix.RData", envir = count_data)
if (!exists("genes_counts", envir = count_data, inherits = FALSE)) {
  stop("The RData file must contain an object named genes_counts.")
}
counts <- count_data$genes_counts

metadata <- read.csv(
  "metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Align metadata with count columns and select cerebellum samples.
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata <- metadata[
  metadata$Sample.Source %in% "Cerebellum",
  ,
  drop = FALSE
]

if (nrow(metadata) == 0L) {
  stop("No cerebellum samples were found in the metadata.")
}
counts <- counts[, rownames(metadata), drop = FALSE]
message("Selected ", nrow(metadata), " cerebellum samples.")

# Match the original conversion of fractional RSEM counts: as.integer()
# truncates fractional values rather than rounding them.
for (i in seq_len(ncol(counts))) {
  counts[, i] <- as.integer(counts[, i])
}

# Remove Ensembl version suffixes and retain the first occurrence of each ID.
counts$gene <- sub("\\..*", "", rownames(counts))
counts <- counts[!duplicated(counts$gene), , drop = FALSE]

# Retrieve HGNC symbols for the bulk Ensembl IDs.
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = counts$gene,
  mart = mart
)

# Retain the first occurrence of each symbol, as in the original analysis.
gene_map <- gene_map[!duplicated(gene_map$hgnc_symbol), , drop = FALSE]
gene_map <- na.omit(gene_map)

# Align bulk counts and gene symbols, removing incomplete rows.
counts <- counts[match(gene_map$ensembl_gene_id, counts$gene), , drop = FALSE]
counts <- counts[!duplicated(counts), , drop = FALSE]
counts <- na.omit(counts)
counts <- counts[counts$gene %in% gene_map$ensembl_gene_id, , drop = FALSE]
gene_map <- gene_map[
  match(counts$gene, gene_map$ensembl_gene_id),
  ,
  drop = FALSE
]
counts$gene <- NULL
rownames(counts) <- gene_map$hgnc_symbol
bulk_matrix <- data.matrix(counts)

# Read the single-cell reference and retain control cells.
reference <- zellkonverter::readH5AD("GSE246183_230502_human_ctype.h5ad")
reference <- reference[, reference$genotype == "Ctrl"]

# Correct reference batch labels, matching the original analysis.
reference$batch[reference$batch == "4092"] <- "409"
reference$batch[reference$batch == "15162"] <- "1516"

# MuSiC uses the reference assay named counts. The original analysis used X.
reference_assays <- SummarizedExperiment::assays(reference)
reference_assays[["counts"]] <- reference_assays[["X"]]
reference_assays[["X"]] <- NULL
SummarizedExperiment::assays(reference) <- reference_assays

saveRDS(reference, "GSE246183_230502_human_ctype_converted.rds")

# Estimate proportions using sub4 cell types and batch as the sample label.
results <- MuSiC::music_prop(
  bulk.mtx = bulk_matrix,
  sc.sce = reference,
  clusters = "sub4",
  samples = "batch",
  verbose = TRUE
)

write.csv(
  results[["Est.prop.weighted"]],
  "music_props_vs_gse246183_only_cerebellum.csv"
)