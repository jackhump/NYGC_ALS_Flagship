# Differential expression analysis: ALS versus control
#
# Run from a working directory containing these files:
#   gene_counts.csv
#   metadata.csv
#   gencode.v30_geneInfo.csv
#
# Install DESeq2 if needed:
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Input requirements:
# - The count CSV contains raw, non-negative integer counts, with gene IDs in
#   the first column and sample IDs as the remaining column names.
# - The metadata CSV has sample IDs in its first column and includes
#   Subject.Group, Sample.Source, Prep, Sex, Age.at.Death, RIN, and
#   Site.Specimen.Collected.
# - The gene annotation CSV has gene IDs in its fifth column and includes
#   biotype and symbol columns.
# - Sample IDs must match between the metadata and count matrix. Gene IDs must
#   match the annotation to obtain gene symbols.
#
# For each tissue, the analysis selects ALS and Control samples, removes
# samples missing age, RIN, or collection site, and fits a DESeq2 model.
# The model adjusts for preparation, sex, age, RIN, and collection site.
# Numeric covariates are not scaled. The reported log2 fold change is
# ALS / Control.
#
# Only lincRNA and protein-coding annotations are used to assign symbols;
# this does not filter genes from the analysis. Unmatched genes have missing
# symbols. Results are written as one CSV per tissue in the working directory.
# Existing output files with the same names will be overwritten.

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. See the installation instructions above.")
}

# Load count data 
counts <- read.csv(
  "gene_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Load metadata
metadata <- read.csv(
  "metadata.csv",
  row.names = 1,
  check.names = FALSE
)
metadata$Disease <- metadata$Subject.Group

# Load gene annotation file
gene_info <- read.csv(
  "gencode.v30_geneInfo.csv",
  header = TRUE,
  row.names = 5
)

# Limit to just protein coding and lincRNA
gene_info <- gene_info[
  gene_info$biotype %in% c("lincRNA", "protein_coding"),
  ,
  drop = FALSE
]
gene_symbols <- setNames(gene_info$symbol, rownames(gene_info))

# Define function to run DEG analysis 
run_deg_analysis <- function(
    tissue,
    count_matrix,
    metadata,
    gene_symbols,
    group_a = "ALS",
    group_b = "Control"
) {
  # Select the tissue and comparison groups.
  metadata <- metadata[
    metadata$Sample.Source == tissue &
      metadata$Disease %in% c(group_a, group_b),
    ,
    drop = FALSE
  ]
  
  # Exclude samples missing these three model covariates.
  metadata <- metadata[
    !is.na(metadata$Age.at.Death) &
      !is.na(metadata$RIN) &
      !is.na(metadata$Site.Specimen.Collected),
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata) == 0L) {
    stop("No eligible samples for tissue: ", tissue)
  }
  
  missing_counts <- setdiff(rownames(metadata), colnames(count_matrix))
  if (length(missing_counts) > 0L) {
    stop(
      "Missing count columns for: ",
      paste(missing_counts, collapse = ", ")
    )
  }
  
  message(tissue, ": ", nrow(metadata), " eligible samples")
  
  design <- ~ Disease + Prep + Sex + Age.at.Death + RIN +
    Site.Specimen.Collected
  
  # Align count columns with the selected metadata rows.
  count_matrix <- count_matrix[, rownames(metadata), drop = FALSE]
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = design
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  results <- data.frame(
    DESeq2::results(dds, contrast = c("Disease", group_a, group_b))
  )
  results$symbol <- unname(gene_symbols[rownames(results)])
  results
}

tissues <- c(
  "Cerebellum",
  "Cortex_Frontal",
  "Spinal_Cord_Cervical",
  "Spinal_Cord_Lumbar"
)

# Run DEG analysis
for (tissue in tissues) {
  results <- run_deg_analysis(
    tissue = tissue,
    count_matrix = counts,
    metadata = metadata,
    gene_symbols = gene_symbols
  )
  
  output_file <- paste0(tissue, "_als_vs_control_simple.csv")
  write.csv(results, output_file, row.names = TRUE)
}

