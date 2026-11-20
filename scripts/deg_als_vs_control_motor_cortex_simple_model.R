# DREAM differential expression analysis: ALS versus control
#
# Run from a working directory containing:
#   metadata_mc.csv
#   genes_counts.csv
#
# Required Bioconductor packages: edgeR, variancePartition, and BiocParallel.
# Install them if needed:
#   install.packages("BiocManager")
#   BiocManager::install(c("edgeR", "variancePartition", "BiocParallel"))
#
# Input requirements:
# - The metadata CSV has sample IDs in its first column and contains Disease,
#   tissue, DNA, Prep, Age.at.Death, Sex, RIN, and Site.Specimen.Collected.
# - The count CSV has gene IDs in its first column and sample IDs as the
#   remaining column names. Counts must be raw, non-negative integer counts.
# - Metadata sample IDs must match count matrix column names.
#
# Analysis:
# - Recode "ALS Spectrum MND" as ALS and every other non-missing Disease value
#   as Control, as in the original analysis.
# - Keep samples with known sex and age, then remove samples missing any model
#   variable. Align count columns to the remaining metadata rows.
# - Calculate edgeR normalization factors and fit a DREAM model with tissue
#   and DNA as random intercepts. The model adjusts for preparation, age,
#   sex, RIN, and collection site.
# - Report the ALS coefficient relative to Control for all genes returned by
#   topTable. No gene filtering is performed in this script.
#
# The output CSV is written to the working directory and contains gene IDs in
# its first column. An existing file with the same name will be overwritten.

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

# Load  metadata
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

metadata <- metadata[
  metadata$Sex %in% c("Female", "Male") &
    !is.na(metadata$Age.at.Death),
  ,
  drop = FALSE
]

# Load counts data
counts <- read.csv(
  "genes_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Check and align samples before selecting complete model records.
missing_counts <- setdiff(rownames(metadata), colnames(counts))
if (length(missing_counts) > 0L) {
  stop(
    "Missing count columns for: ",
    paste(missing_counts, collapse = ", ")
  )
}
counts <- counts[, rownames(metadata), drop = FALSE]

# Select model variables and remove samples missing any of them.
model_variables <- c(
  "Disease", "tissue", "DNA", "Prep", "Age.at.Death",
  "Sex", "RIN", "Site.Specimen.Collected"
)
model_metadata <- metadata[, model_variables, drop = FALSE]
model_metadata <- na.omit(model_metadata)

if (nrow(model_metadata) == 0L) {
  stop("No samples remain after filtering for complete model variables.")
}

counts <- counts[, rownames(model_metadata), drop = FALSE]
model_metadata$Disease <- factor(
  model_metadata$Disease,
  levels = c("Control", "ALS")
)

message(
  "Model data: ", nrow(model_metadata), " samples and ",
  nrow(counts), " genes"
)

# Calculate normalization factors for the count data.
dge <- edgeR::DGEList(counts)
dge <- edgeR::calcNormFactors(dge)

# Tissue and DNA account for grouping through random intercepts.
model_formula <- ~ Disease + Prep + Age.at.Death + Sex + RIN +
  Site.Specimen.Collected + (1 | tissue) + (1 | DNA)

# Use four socket workers, matching the original analysis.
parallel_parameters <- BiocParallel::SnowParam(
  workers = 4,
  type = "SOCK",
  progressbar = TRUE
)
BiocParallel::register(parallel_parameters)

# Calculate expression weights and fit the mixed model.
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

# Extract the ALS-versus-Control coefficient for all available genes.
results <- variancePartition::topTable(
  fit,
  coef = "DiseaseALS",
  number = Inf
)

write.csv(
  results,
  "Cortex_Motor_als_vs_control_simple.csv",
  row.names = TRUE
)

