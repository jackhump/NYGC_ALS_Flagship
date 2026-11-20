# DESeq2 differential expression analysis: SOD1-positive versus SOD1-negative
#
# Run from a working directory containing:
#   metadata.csv
#   gene_counts.csv
#   gencode.v30_geneInfo.csv
#
# Required packages: readr and DESeq2. Install them if needed:
#   install.packages("readr")
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Input requirements:
# - The covariate TSV has sample IDs in X and includes Sample.Source, sod1,
#   Sex, Age.at.Symptom.Onset, and Age.at.Death.
# - The count CSV has gene IDs in its first column and sample IDs as the
#   remaining column names. Counts must be raw, non-negative integers.
# - The annotation CSV has gene IDs in its fifth column and includes biotype
#   and symbol columns. Gene IDs must match the count data to obtain symbols.
# - Metadata sample IDs must match count matrix column names.
#
# For each tissue, the script selects samples with known ages and sod1 values
# of Yes or No. It fits a DESeq2 model adjusted for sex, age at symptom onset,
# and age at death. The reported log2 fold change is Yes / No.
#
# Annotation is restricted to lincRNA and protein-coding symbols; this does not
# filter genes from the analysis. Unmatched genes have missing symbols.
# The function also identifies genes with adjusted p-value below 0.05, but the
# output files contain the full results, including genes with missing p-values.
# One CSV per tissue is written to the working directory, with gene IDs in the
# first column. Existing output files with the same names will be overwritten.

required_packages <- c("readr", "DESeq2")
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

metadata <- as.data.frame(readr::read_tsv("metadata.csv"))
rownames(metadata) <- metadata$X

# Load gene counts data
counts <- read.csv(
  "gene_counts.csv",
  row.names = 1,
  check.names = FALSE
)

# Exclude samples missing either age measure.
metadata <- metadata[
  !is.na(metadata$Age.at.Death) &
    !is.na(metadata$Age.at.Symptom.Onset),
  ,
  drop = FALSE
]

# Load gene annotation file 
gene_info <- read.csv(
  "gencode.v30_geneInfo.csv",
  header = TRUE,
  row.names = 5
)

# Limit to just protein_coding and lincRNA
gene_info <- gene_info[
  gene_info$biotype %in% c("lincRNA", "protein_coding"),
  ,
  drop = FALSE
]
gene_symbols <- setNames(gene_info$symbol, rownames(gene_info))

# Define DEG function
run_deg_analysis <- function(
    tissue,
    count_matrix,
    metadata,
    gene_symbols,
    target,
    group_a,
    group_b,
    threshold = 0.05
) {
  # Select the tissue and the two groups being compared.
  sample_metadata <- metadata[
    metadata$Sample.Source == tissue &
      metadata[[target]] %in% c(group_a, group_b),
    ,
    drop = FALSE
  ]
  
  if (nrow(sample_metadata) == 0L) {
    stop("No eligible samples for tissue: ", tissue)
  }
  
  missing_counts <- setdiff(rownames(sample_metadata), colnames(count_matrix))
  if (length(missing_counts) > 0L) {
    stop(
      "Missing count columns for: ",
      paste(missing_counts, collapse = ", ")
    )
  }
  count_matrix <- count_matrix[, rownames(sample_metadata), drop = FALSE]
  
  # Adjust for sex and both age measures, as in the original analysis.
  design <- reformulate(
    c(target, "Sex", "Age.at.Symptom.Onset", "Age.at.Death")
  )
  
  message(tissue, " group counts:")
  print(table(sample_metadata[[target]]))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = design
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  results <- data.frame(
    DESeq2::results(dds, contrast = c(target, group_a, group_b))
  )
  results$symbol <- unname(gene_symbols[rownames(results)])
  
  significant_results <- results[
    !is.na(results$padj) & results$padj < threshold,
    ,
    drop = FALSE
  ]
  
  list(results_full = results, results_DE = significant_results)
}

tissues <- c(
  "Cerebellum",
  "Cortex_Frontal",
  "Spinal_Cord_Cervical",
  "Spinal_Cord_Lumbar"
)

# Run DEG per tissue 
for (tissue in tissues) {
  deg_results <- run_deg_analysis(
    tissue = tissue,
    count_matrix = counts,
    metadata = metadata,
    gene_symbols = gene_symbols,
    target = "sod1",
    group_a = "Yes",
    group_b = "No"
  )
  
  output_file <- paste0(tissue, "_sod1_pos_vs_neg.csv")
  write.csv(deg_results$results_full, output_file, row.names = TRUE)
  
}

