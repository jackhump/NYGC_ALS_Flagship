# Estimate cortex cell-type proportions with MuSiC
#
# Run from a working directory containing:
#   gene_matrix.RData
#   metadata.csv
#   genemes_names_df.csv
#   mathys_all_cogdx_control_cortex_scref_seurat_object.rds
#
# Required packages: Seurat and MuSiC. The original environment used R 4.3.1
# and Seurat 4.0.1; use compatible package versions and dependencies.
#
# Input requirements:
# - The RData file contains genes_counts, a data frame with versioned gene IDs
#   as row names and bulk RNA sample IDs as column names.
# - The metadata CSV has RNA sample IDs in its first column and a Sample.Source
#   column. All count samples are assumed to have matching metadata.
# - The gene-name CSV contains gene_ids and genes columns.
# - The RDS file contains a Seurat reference with broad.cell.type and projid
#   fields in its cell metadata.
#
# The script selects motor and frontal cortex samples, matches bulk genes to
# reference gene names, and estimates cell-type proportions with MuSiC.
# The weighted proportion matrix is written to the working directory as a CSV.
# An existing output file with the same name will be overwritten.

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

# Align metadata with the count matrix and select the two cortex tissues.
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata <- metadata[
  metadata$Sample.Source %in% c("Cortex_Motor", "Cortex_Frontal"),
  ,
  drop = FALSE
]

if (nrow(metadata) == 0L) {
  stop("No motor or frontal cortex samples were found in the metadata.")
}
counts <- counts[, rownames(metadata), drop = FALSE]
message("Selected ", nrow(metadata), " cortex samples.")

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

# Convert the cortex reference to a SingleCellExperiment and run MuSiC.
single_cell_reference <- readRDS(
  "mathys_all_cogdx_control_cortex_scref_seurat_object.rds"
)
single_cell_reference <- Seurat::as.SingleCellExperiment(single_cell_reference)

results <- MuSiC::music_prop(
  bulk.mtx = bulk_matrix,
  sc.sce = single_cell_reference,
  clusters = "broad.cell.type",
  samples = "projid",
  verbose = TRUE
)

write.csv(
  results[["Est.prop.weighted"]],
  "music_props_vs_mathys_motor_and_frontal_cortex.csv"
)

