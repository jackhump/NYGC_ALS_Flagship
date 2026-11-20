# DREAM differential expression analysis: ALS versus control
#
# Run from a working directory containing:
#   metadata_mc.csv
#   gene_counts.csv
#
# Required Bioconductor packages: edgeR, variancePartition, and BiocParallel.
# Install them if needed:
#   install.packages("BiocManager")
#   BiocManager::install(c("edgeR", "variancePartition", "BiocParallel"))
#
# Input requirements:
# - The metadata CSV has sample IDs in its first column and includes Disease,
#   MEAN_READ_LENGTH, Prep, PCT_CODING_BASES, DNA, and tissue.
# - The count CSV has gene IDs in its first column and sample IDs as the
#   remaining column names. Counts must be raw, non-negative integers.
# - Every metadata sample ID must match a count matrix column name.
#
# Analysis:
# - Recode "ALS Spectrum MND" as ALS and every other non-missing Disease value
#   as Control. Missing disease values remain missing.
# - Standardize mean read length and percentage of coding bases across samples.
# - Calculate edgeR normalization factors and fit a DREAM model adjusted for
#   mean read length, preparation, and percentage of coding bases, with tissue
#   and DNA random intercepts.
# - Report the ALS coefficient relative to Control for all genes returned by
#   topTable. No gene-expression filter is applied.
#
# The output CSV is written to the working directory with gene IDs in its first
# column. An existing file with the same name will be overwritten.

required_packages <- c("edgeR", "variancePartition", "BiocParallel")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". See the installation instructions above."
  )
}

metadata <- read.csv(
  "metadata_mc.csv",
  row.names = 1,
  check.names = FALSE
)
metadata$Disease <- ifelse(
  metadata$Disease == "ALS Spectrum MND",
  "ALS",
  "Control"
)

counts <- read.csv(
  "gene_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Match count columns to metadata sample order.
missing_counts <- setdiff(rownames(metadata), colnames(counts))
if (length(missing_counts) > 0L) {
  stop(
    "Missing count columns for: ",
    paste(missing_counts, collapse = ", ")
  )
}
counts <- counts[, rownames(metadata), drop = FALSE]

model_metadata <- metadata[
  ,
  c("Disease", "MEAN_READ_LENGTH", "Prep", "PCT_CODING_BASES", "DNA", "tissue"),
  drop = FALSE
]
model_metadata$Disease <- factor(
  model_metadata$Disease,
  levels = c("Control", "ALS")
)

# Standardize the two numeric covariates.
columns_to_scale <- c("MEAN_READ_LENGTH", "PCT_CODING_BASES")
model_metadata[columns_to_scale] <- lapply(
  model_metadata[columns_to_scale],
  function(x) as.numeric(scale(x))
)

# Calculate normalization factors using all genes.
dge <- edgeR::DGEList(counts)
dge <- edgeR::calcNormFactors(dge)

model_formula <- ~ Disease + MEAN_READ_LENGTH + Prep + PCT_CODING_BASES +
  (1 | tissue) + (1 | DNA)

# Use four socket workers, matching the original analysis.
parallel_parameters <- BiocParallel::SnowParam(
  workers = 4,
  type = "SOCK",
  progressbar = TRUE
)
BiocParallel::register(parallel_parameters)

weighted_expression <- variancePartition::voomWithDreamWeights(
  dge,
  model_formula,
  model_metadata
)
fit <- variancePartition::dream(
  weighted_expression,
  model_formula,
  model_metadata
)

# Extract the ALS-versus-Control coefficient.
results <- variancePartition::topTable(
  fit,
  coef = "DiseaseALS",
  number = Inf
)

write.csv(
  results,
  "Cortex_Motor_als_vs_ctrl.csv",
  row.names = TRUE
)

