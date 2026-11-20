# Estimate spinal cord cell-type proportions with MuSiC
#
# Run from a working directory containing:
#   gene_matrix.RData
#   metadata.csv
#   genemes_names_df.csv
#   human_spinalcord_GSE190442.rds
#
# Required packages: Seurat and MuSiC. The original environment used R 4.3.1
# and Seurat 4.0.1; use compatible package versions and dependencies.
#
# Input requirements:
# - The RData file contains genes_counts, a data frame with versioned gene IDs
#   as row names and bulk RNA sample IDs as column names.
# - metadata.csv has RNA sample IDs in its first column and a Sample.Source
#   column identifying the tissue. All count samples have matching metadata.
# - genemes_names_df.csv contains gene_ids and genes columns.
# - The RDS file contains a Seurat reference with top_level_annotation and
#   sample fields in its cell metadata.
#
# The script selects all lumbar, cervical, thoracic, and unspecified spinal
# cord samples, matches bulk genes to reference gene names, and estimates
# cell-type proportions with MuSiC. DNA matching is not required.
#
# The complete MuSiC result and its weighted proportion matrix are written to
# the working directory as RDS and CSV files. Existing files are overwritten.

required_packages <- c("Seurat", "MuSiC")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

# Load counts separately so unrelated objects in the RData file do not enter
# the main workspace.
count_data <- new.env()
load("gene_matrix.RData", envir = count_data)
if (!exists("genes_counts", envir = count_data, inherits = FALSE)) {
  stop("The RData file must contain an object named genes_counts.")
}
counts <- count_data$genes_counts

metadata <- read.csv(
  "metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Align metadata rows with count columns.
metadata <- metadata[colnames(counts), , drop = FALSE]

spinal_cord_tissues <- c(
  "Spinal_Cord_Lumbar",
  "Spinal_Cord_Cervical",
  "Spinal_Cord_Thoracic",
  "Spinal_Cord_Unspecified"
)

# Trim whitespace to recognize the unspecified label with a trailing space.
metadata <- metadata[
  trimws(metadata$Sample.Source) %in% spinal_cord_tissues,
  ,
  drop = FALSE
]

if (nrow(metadata) == 0L) {
  stop("No spinal cord samples were found in the metadata.")
}
counts <- counts[, rownames(metadata), drop = FALSE]
message("Selected ", nrow(metadata), " spinal cord samples.")

# Match the original conversion of fractional RSEM counts: as.integer()
# truncates fractional values rather than rounding them.
for (i in seq_len(ncol(counts))) {
  counts[, i] <- as.integer(counts[, i])
}

# Remove Ensembl version suffixes and retain the first occurrence of each ID.
counts$gene <- sub("\\..*$", "", rownames(counts))
counts <- counts[!duplicated(counts$gene), , drop = FALSE]

gene_names <- read.csv("genemes_names_df.csv")
gene_names$X <- NULL
gene_names <- na.omit(gene_names)

# Match bulk genes to the gene-name table and remove incomplete rows.
counts <- counts[match(gene_names$gene_ids, counts$gene), , drop = FALSE]
counts <- na.omit(counts)
gene_names <- gene_names[
  match(counts$gene, gene_names$gene_ids),
  ,
  drop = FALSE
]
counts$gene <- NULL
rownames(counts) <- gene_names$genes
bulk_matrix <- data.matrix(counts)

# Convert the reference to a SingleCellExperiment and run MuSiC.
single_cell_reference <- readRDS("human_spinalcord_GSE190442.rds")
single_cell_reference <- Seurat::as.SingleCellExperiment(single_cell_reference)

results <- MuSiC::music_prop(
  bulk.mtx = bulk_matrix,
  sc.sce = single_cell_reference,
  clusters = "top_level_annotation",
  samples = "sample",
  verbose = TRUE
)

saveRDS(
  results,
  "music_props_vs_spinal_cord_gse190442_spinal_cord_samples.rds"
)
write.csv(
  results[["Est.prop.weighted"]],
  "music_props_vs_spinal_cord_gse190442_spinal_cord_samples.csv"
)

