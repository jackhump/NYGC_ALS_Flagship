# DREAM differential expression analysis: SOD1-positive versus SOD1-negative
#
# Run from a working directory containing:
#   covs_with_sod1.tsv
#   metadata_mc.csv
#   genes_counts.csv
#
# Required Bioconductor packages: edgeR, variancePartition, and BiocParallel.
# Install them if needed:
#   install.packages("BiocManager")
#   BiocManager::install(c("edgeR", "variancePartition", "BiocParallel"))
#
# Input requirements:
# - covs_with_sod1.tsv has sample IDs in X and includes Sample.Source, sod1,
#   Prep, DNA, Sex, Age.at.Death, Age.at.Symptom.Onset,
#   Disease.Duration.in.Months, and Site.of.Motor.Onset.
# - metadata_mc.csv has sample IDs in its first column and a
#   tissue column. It supplies the tissue grouping variable for the model.
# - genes_counts.csv has gene IDs in its first column and sample IDs as
#   the remaining column names. Counts must be raw, non-negative integers.
# - Sample IDs must match across the files. Age and duration variables must be
#   numeric, with missing values represented as NA.
#
# Analysis:
# - Attach the tissue variable and select Cortex_Motor samples.
# - When possible, derive missing disease duration, onset age, or age at death
#   from the other two values. Keep samples with both ages available.
# - Calculate edgeR normalization factors and fit a DREAM model comparing
#   sod1 = Yes with sod1 = No. The model adjusts for preparation, age at death,
#   age at symptom onset, and sex, with tissue and DNA random intercepts.
# - No gene-expression filter is applied; all genes are passed to the model.
#
# The output CSV is written to the working directory, with gene IDs in its
# first column. An existing output file with the same name will be overwritten.

required_packages <- c("readr", "edgeR", "variancePartition", "BiocParallel")
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

# Load covariates and attach the original tissue grouping.
metadata <- as.data.frame(readr::read_tsv("covs_with_sod1.tsv"))
metadata_orig <- read.csv(
  "metadata_mc.csv",
  row.names = 1,
  check.names = FALSE
)
metadata_orig$sample <- rownames(metadata_orig)
metadata_orig <- metadata_orig[, c("sample", "tissue"), drop = FALSE]

metadata <- merge(
  metadata,
  metadata_orig,
  by.x = "X",
  by.y = "sample",
  all.x = TRUE
)
rownames(metadata) <- metadata$X
metadata <- metadata[
  metadata$Sample.Source == "Cortex_Motor",
  ,
  drop = FALSE
]

# Fill missing values when the other two related measurements are available.
missing_duration <- is.na(metadata$Disease.Duration.in.Months)
metadata$Disease.Duration.in.Months[missing_duration] <-
  (metadata$Age.at.Death[missing_duration] -
     metadata$Age.at.Symptom.Onset[missing_duration]) * 12

missing_onset <- is.na(metadata$Age.at.Symptom.Onset)
metadata$Age.at.Symptom.Onset[missing_onset] <-
  metadata$Age.at.Death[missing_onset] -
  metadata$Disease.Duration.in.Months[missing_onset] / 12

missing_death <- is.na(metadata$Age.at.Death)
metadata$Age.at.Death[missing_death] <-
  metadata$Age.at.Symptom.Onset[missing_death] +
  metadata$Disease.Duration.in.Months[missing_death] / 12

metadata <- metadata[
  !is.na(metadata$Age.at.Death) &
    !is.na(metadata$Age.at.Symptom.Onset),
  ,
  drop = FALSE
]

if (nrow(metadata) == 0L) {
  stop("No motor cortex samples remain after filtering for age.")
}

counts <- read.csv(
  "genes_counts.csv",
  row.names = 1,
  check.names = FALSE
)
missing_counts <- setdiff(rownames(metadata), colnames(counts))
if (length(missing_counts) > 0L) {
  stop(
    "Missing count columns for: ",
    paste(missing_counts, collapse = ", ")
  )
}
counts <- counts[, rownames(metadata), drop = FALSE]

# Preserve the original set of metadata columns supplied to DREAM.
model_metadata <- metadata[
  ,
  c(
    "Prep", "DNA", "tissue", "Age.at.Death", "Age.at.Symptom.Onset",
    "Site.of.Motor.Onset", "Sex", "Disease.Duration.in.Months",
    "sod1", "Sample.Source"
  ),
  drop = FALSE
]
model_metadata$sod1 <- factor(model_metadata$sod1, levels = c("No", "Yes"))

message("Motor cortex SOD1 sample counts:")
print(table(model_metadata$sod1, useNA = "ifany"))

# Calculate normalization factors using all genes.
dge <- edgeR::DGEList(counts)
dge <- edgeR::calcNormFactors(dge)

model_formula <- ~ sod1 + Prep + Age.at.Death + Age.at.Symptom.Onset +
  Sex + (1 | tissue) + (1 | DNA)

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

# Extract the SOD1-positive versus SOD1-negative coefficient.
results <- variancePartition::topTable(
  fit,
  coef = "sod1Yes",
  number = Inf
)

write.csv(
  results,
  "Cortex_Motor_sod1_pos_vs_neg.csv",
  row.names = TRUE
)

