# Differential expression analysis: ALS versus control, adjusted for cell types
#
# Run from the directory containing the input files:
#   Rscript cleaned_deg_analysis.R
# Install the Bioconductor package DESeq2 before running this script, for example:
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Required files in the working directory:
#   gene_count.csv: raw gene counts; first column is gene ID,
#     remaining columns are sample IDs.
#   metadata.csv: sample metadata; first column is sample ID.
#   props_cerebellum_norm.tsv, props_cortex_norm.tsv, and
#     props_spinal_cord_norm.tsv: normalized cell-type proportions, with a
#     Sample column containing sample IDs and remaining columns as covariates.
#   gencode.v30_geneInfo.csv: gene annotation; fifth column is gene ID, and
#     columns named biotype and symbol provide gene type and gene symbol.
#
# Sample IDs must match across the count matrix, metadata, and the applicable
# cell-type file. Gene IDs must match the annotation IDs to obtain symbols.
# The count matrix must contain non-negative integer raw counts suitable for
# DESeq2. Metadata must include Subject.Group, Sample.Source, and the covariates
# used in the tissue-specific designs below. Each design must be estimable.
#
# For each tissue, this script selects ALS and Control samples, standardizes
# numeric metadata columns within that subset, adds the corresponding cell-type
# covariates, and fits DESeq2. The reported log2 fold change is ALS / Control.
# Annotation is restricted to lincRNA and protein_coding symbols; this does not
# filter genes from the differential expression analysis. Unmatched genes have
# missing symbols. One CSV per tissue is written to the working directory, with
# gene IDs in the first column. Existing output files will be overwritten.

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("The DESeq2 package is required; see the installation instructions above.")
}

# Load gene counts data
counts <- read.csv(
  "gene_count.csv", row.names = 1, check.names = FALSE
)

# Load metadata
metadata <- read.csv(
  "metadata.csv", row.names = 1, check.names = FALSE
)
metadata$Disease <- metadata$Subject.Group

# Load cell type proportions
# Different files as single-cell reference differed per tissue source
cell_type_files <- c(
  cerebellum = "props_cerebellum_norm.tsv",
  cortex = "props_cortex_norm.tsv",
  spinal_cord = "props_spinal_cord_norm.tsv"
)

# Convert to z-scores
cell_type_z_scores <- lapply(cell_type_files, function(file) {
  read.delim(file, row.names = "Sample", check.names = FALSE)
})

# Both spinal cord regions use the same cell-type proportion file.
cell_type_source <- c(
  Cerebellum = "cerebellum",
  Cortex_Frontal = "cortex",
  Spinal_Cord_Cervical = "spinal_cord",
  Spinal_Cord_Lumbar = "spinal_cord"
)

# Load gene annotation file
gene_info <- read.csv("gencode.v30_geneInfo.csv", header = TRUE, row.names = 5)

# Limit to just protein-coding and lincRNA
gene_info <- gene_info[gene_info$biotype %in% c("lincRNA", "protein_coding"), ]

# Map ensembl IDs to gene symbols 
gene_symbols <- setNames(gene_info$symbol, rownames(gene_info))

# Define function to run differential gene analysis 
run_deg_analysis <- function(tissue, count_matrix, metadata, design,
                             cell_type_scores, gene_symbols,
                             group_a = "ALS", group_b = "Control") {
  # Restrict samples before scaling so each tissue is standardized separately.
  metadata <- metadata[
    metadata$Sample.Source == tissue & metadata$Disease %in% c(group_a, group_b),
    , drop = FALSE
  ]
  if (nrow(metadata) == 0L) {
    stop("No matching samples for tissue: ", tissue)
  }
  numeric_columns <- vapply(metadata, is.numeric, logical(1))
  metadata[numeric_columns] <- lapply(metadata[numeric_columns], function(x) {
    as.numeric(scale(x))
  })
  
  # Match by sample ID rather than assuming the files share the same order.
  score_rows <- match(rownames(metadata), rownames(cell_type_scores))
  if (anyNA(score_rows)) {
    stop("Missing cell-type scores for: ",
         paste(rownames(metadata)[is.na(score_rows)], collapse = ", "))
  }
  missing_counts <- setdiff(rownames(metadata), colnames(count_matrix))
  if (length(missing_counts) > 0L) {
    stop("Missing count columns for: ", paste(missing_counts, collapse = ", "))
  }
  cell_type_scores <- cell_type_scores[score_rows, , drop = FALSE]
  rownames(cell_type_scores) <- rownames(metadata)
  metadata <- cbind(metadata, cell_type_scores)
  
  # Add every cell-type column to the tissue-specific adjustment formula.
  design <- reformulate(
    c(attr(terms(design), "term.labels"), colnames(cell_type_scores)),
    env = environment(design)
  )
  count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix, colData = metadata, design = design
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  results <- data.frame(DESeq2::results(dds, contrast = c("Disease", group_a, group_b)))
  results$symbol <- unname(gene_symbols[rownames(results)])
  results
}

# Covariates differ by tissue, matching the original analysis.
designs <- list(
  Cerebellum = ~ Disease + Prep + MEAN_READ_LENGTH + Site.Specimen.Collected,
  Cortex_Frontal = ~ Disease + Prep + MEAN_READ_LENGTH + PCT_CHIMERAS,
  Spinal_Cord_Cervical = ~ Disease + Prep + MEAN_READ_LENGTH + Site.Specimen.Collected,
  Spinal_Cord_Lumbar = ~ Disease + Site.Specimen.Collected + Prep + PCT_CHIMERAS
)


# Run differential gene expression per tissue 
for (tissue in names(designs)) {
  results <- run_deg_analysis(
    tissue = tissue,
    count_matrix = counts,
    metadata = metadata,
    design = designs[[tissue]],
    cell_type_scores = cell_type_z_scores[[cell_type_source[[tissue]]]],
    gene_symbols = gene_symbols
  )
  output_file <- paste0(tissue, "_als_vs_control_plus_cell_types.csv")
  write.csv(results, output_file)
}
