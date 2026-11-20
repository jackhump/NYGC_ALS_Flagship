# DESeq2 differential expression analysis: ALS versus control
#
# Run from a working directory containing:
#   gene_counts.csv
#   metadata.csv
#   gencode.v30_geneInfo.csv
#
# Install DESeq2 if needed:
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Input requirements:
# - The count CSV has gene IDs in its first column and sample IDs as the
#   remaining column names. Counts must be raw, non-negative integers.
# - The metadata CSV has sample IDs in its first column and includes
#   Subject.Group, Sample.Source, and the covariates in the selected design.
# - The annotation CSV has gene IDs in its fifth column and includes biotype
#   and symbol columns.
# - Metadata sample IDs must match count matrix column names. Gene IDs must
#   match the annotation to obtain symbols.
#
# Analysis:
# - Select ALS and Control samples for each tissue listed in `tissues`.
# - Standardize numeric metadata columns within each tissue subset.
# - Fit the tissue-specific DESeq2 model and report ALS / Control results.
# - No gene-expression filter is applied. Annotation is restricted to lincRNA
#   and protein-coding symbols, but this does not filter genes from the analysis.
# - The function identifies genes with adjusted p-value below 0.05; the CSV
#   contains the full results, including genes with missing adjusted p-values.
#
# Results are written to the working directory with gene IDs in the first
# column. Existing output files with the same names will be overwritten.

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required. See the installation instructions above.")
}

counts <- read.csv(
  "gene_counts.csv",
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  "metadata.csv",
  row.names = 1,
  check.names = FALSE
)
metadata$Disease <- metadata$Subject.Group

gene_info <- read.csv(
  "gencode.v30_geneInfo.csv",
  header = TRUE,
  row.names = 5
)
gene_info <- gene_info[
  gene_info$biotype %in% c("lincRNA", "protein_coding"),
  ,
  drop = FALSE
]
gene_symbols <- setNames(gene_info$symbol, rownames(gene_info))

run_deg_analysis <- function(
    tissue,
    count_matrix,
    metadata,
    design,
    gene_symbols,
    threshold = 0.05,
    group_a = "ALS",
    group_b = "Control"
) {
  # Select the tissue and comparison groups.
  sample_metadata <- metadata[
    metadata$Sample.Source == tissue &
      metadata$Disease %in% c(group_a, group_b),
    ,
    drop = FALSE
  ]
  
  if (nrow(sample_metadata) == 0L) {
    stop("No eligible samples for tissue: ", tissue)
  }
  
  # Standardize numeric metadata within this tissue subset.
  numeric_columns <- vapply(sample_metadata, is.numeric, logical(1))
  sample_metadata[numeric_columns] <- lapply(
    sample_metadata[numeric_columns],
    function(x) as.numeric(scale(x))
  )
  
  missing_counts <- setdiff(rownames(sample_metadata), colnames(count_matrix))
  if (length(missing_counts) > 0L) {
    stop(
      "Missing count columns for: ",
      paste(missing_counts, collapse = ", ")
    )
  }
  count_matrix <- count_matrix[, rownames(sample_metadata), drop = FALSE]
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = design
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  
  results <- data.frame(
    DESeq2::results(dds, contrast = c("Disease", group_a, group_b))
  )
  results$symbol <- unname(gene_symbols[rownames(results)])
  
  significant_results <- results[
    !is.na(results$padj) & results$padj < threshold,
    ,
    drop = FALSE
  ]
  
  list(
    results_full = results,
    results_DE = significant_results,
    genes_analyzed = nrow(count_matrix)
  )
}


# Define tissues 
tissues <- c("Cerebellum", "Cortex_Frontal", "Cortex_Temporal",
             "Hippocampus", "Spinal_Cord_Cervical", "Spinal_Cord_Lumbar",
             "Spinal_Cord_Thoracic")

designs <- list(
  Cerebellum = ~ Disease + Prep + MEAN_READ_LENGTH + Site.Specimen.Collected,
  Spinal_Cord_Lumbar = ~ Disease + Site.Specimen.Collected + Prep + PCT_CHIMERAS,
  Spinal_Cord_Cervical = ~ Disease + Prep + MEAN_READ_LENGTH + Site.Specimen.Collected,
  Spinal_Cord_Thoracic = ~ Disease + Prep + MEDIAN_5PRIME_BIAS + PCT_CHIMERAS,
  Cortex_Frontal = ~ Disease + Prep + MEAN_READ_LENGTH + PCT_CHIMERAS,
  Hippocampus = ~ Disease + MEDIAN_5PRIME_BIAS + MEDIAN_3PRIME_BIAS + PCT_CODING_BASES,
  Cortex_Motor = ~ Disease + Prep + MEAN_READ_LENGTH + PCT_CHIMERAS,
  Cortex_Temporal = ~ Disease + Prep + PCT_CHIMERAS + MEDIAN_5PRIME_BIAS
)

for (tissue in tissues) {
  if (is.null(designs[[tissue]])) {
    stop("No design defined for tissue: ", tissue)
  }
  
  deg_results <- run_deg_analysis(
    tissue = tissue,
    count_matrix = counts,
    metadata = metadata,
    design = designs[[tissue]],
    gene_symbols = gene_symbols
  )
  
  output_file <- paste0(tissue, "_als_vs_control.csv")
  write.csv(deg_results$results_full, output_file, row.names = TRUE)
}
