# VCF QC Annotation and Filtering Pipeline

This document describes the workflow used to annotate a cohort-level WGS VCF with quality control (QC) metrics and to derive a high-confidence subset of variants based on stringent QC filtering criteria.

---

## Overview

The workflow consists of the following steps:

1. Harmonization and merging of variant callsets (SNV/INDEL and SV)
2. Variant normalization and annotation with QC metrics (Hardy–Weinberg equilibrium and low-complexity regions)
3. Derivation of variant-level QC flags
4. Annotation of the VCF with QC flags and extraction of a high-confidence variant subset

---

## 1. Harmonization and Merging of Variant Callsets

Structural variant (SV) callsets were first harmonized to match the sample ordering of SNV/INDEL data. This was performed on a per-chromosome basis to ensure consistency across variant types.

Sample reordering and chromosome subsetting were performed using:

```bash
bcftools view \
  -r ${chr} \
  --samples-file ${samples_order_list} \
  --threads 10 \
  -Oz -o ${vcf_prefix}.${chr}.samples_reordered.vcf.gz \
  ${vcf}
```

The reordered SV callset was indexed using tabix.

SNV/INDEL and SV callsets were subsequently concatenated and sorted to produce a unified cohort-level VCF:

```bash
bcftools concat \
  --allow-overlaps \
  -Oz -o ${out_prefix}.unsorted.vcf.gz \
  --threads 16 \
  ${vcf1} ${vcf2}

bcftools sort \
  -T ${tmpdir} \
  -m 2G \
  -Oz -o ${out_prefix}.vcf.gz \
  ${out_prefix}.unsorted.vcf.gz
```

## 2. Annotation of QC Metrics & normalization

### 2.1 Hardy–Weinberg Equilibrium

Hardy–Weinberg equilibrium (HWE) p-values were computed using the bcftools +fill-tags plugin:

```bash
bcftools +fill-tags \
  ${vcf} \
  --threads 10 \
  -Oz -o ${prefix}.fill-tags.vcf.gz \
  -- -t HWE
```
### 2.2 Variant normalization

Variants were left-normalized relative to the reference genome and multi-allelic sites were decomposed into bi-allelic records to ensure consistent representation:

```bash
bcftools norm \
  --fasta-ref ${ref_fasta} \
  --multiallelics -both \
  -Oz -o ${prefix}.norm.vcf.gz \
  --threads 10 \
  ${vcf}
```

### 2.3 Low-Complexity Regions

Variants overlapping low-complexity regions were identified using Genome in a Bottle (GIAB) stratification bed files. Variants were flagged if at least 80% of their length overlapped a low-complexity region.

annotation_bed="GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz"; downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/LowComplexity/)

```bash
bcftools annotate \
--annotations ${annotation_bed} \
--mark-sites +LCR80 \
--columns "CHROM,FROM,TO" \
--min-overlap :0.8 \
--threads 10 \
--header-lines ${header_new_lines} \
-Oz -o ${outdir}/${prefix}.LCR_minOverlap80.vcf.gz \
${vcf}
```

## 3. Derivation of QC Flags

Variant-level QC metrics were extracted and combined to define filtering flags. The following metrics were used:

- **Call rate (CR):** proportion of samples with non-missing genotype calls  
- **Hardy–Weinberg equilibrium (HWE) p-value**  
- **Low-complexity region overlap (LCR80)**  

Variants were classified according to the following criteria:

- **Call rate > 0.95**  
- **HWE p-value > 1 × 10⁻¹⁰**  
- **No overlap with low-complexity regions (LCR80 flag absent)**  

These criteria were combined to define a **Stringent_QC** flag.
Conceptually, the filtering logic can be summarized as:

```text
IF call_rate > 0.95 → CR95_QC = PASS ELSE FAIL
IF HWE > 1e-10 → HWE_QC = PASS ELSE FAIL
IF LCR80 present → LCR80_QC = FAIL ELSE PASS

IF CR95_QC == PASS AND HWE_QC == PASS AND LCR80_QC == PASS
    → Stringent_QC = PASS
ELSE
    → Stringent_QC = FAIL
```

In practice, these flags were derived using bcftools query to extract relevant fields, followed by post-processing with awk.

Chromosome-Specific Handling

For chromosome Y, call rate was inferred from allele number (AN), normalized by the total number of male samples:

```text
call_rate = AN / number_of_male_samples
```

This accounts for haploid genotype representation in male samples.

## 4. Annotation with QC Flags & filtering

### 4.1 Annotation of the VCF with QC flags

The derived QC flags (CR95_QC, HWE_QC, LCR80_QC, and Stringent_QC) were added back to the VCF using bcftools annotate:

```bash
bcftools annotate \
--annotations ${annotation_file} \
--columns "CHROM,POS,~ID,REF,ALT,CR95_QC,HWE_QC,LCR80_QC,Stringent_QC" \
--threads 10 \
--pair-logic exact \
--header-lines ${header_new_lines} \
-Oz \
-o ${outdir}/${prefix}.flagged.vcf.gz \
${vcf}
```

### 4.2. Final cleanup & filtering

#### 4.2.1 Removal of Spanning Deletion Alleles

Variants containing spanning deletion alleles (“*”) were removed since they are meaningless after multiallelics splitting:

```bash
bcftools view \
--threads 10 \
--no-version \
 ${vcf} \
| awk '{if($5 != "*") print $0}' | bgzip -c > ${outdir}/${prefix}.clean.vcf.gz
```

#### 4.2.2 Extraction of High-Confidence Variants

Variants passing all QC criteria were extracted to generate the high-confidence callset:

```bash
bcftools view \
--threads 10 \
-i 'INFO/Stringent_QC == "PASS"' \
--write-index \
-Oz \
-o ${outdir}/${prefix}.StringentQC.vcf.gz \
${vcf}
```

# Reproducibility Notes

- The workflow was executed using bcftools (v1.19 or higher) and standard UNIX utilities.
- Execution was performed in a high-performance computing environment; some steps assume parallelization and intermediate file handling.

# Data Dependencies

- Reference genome: GRCh38 (ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
- Low-complexity region annotations: GRCh38 Genome in a Bottle (GIAB) genome stratification bed files (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/LowComplexity/)

# Citations

1. Danecek P, Bonfield JK, et al. Twelve years of SAMtools and BCFtools. Gigascience (2021) 10(2):giab008
2. Dwarshuis, N., Kalra, D., McDaniel, J. et al. The GIAB genomic stratifications resource for human reference genomes. Nat Commun 15, 9029 (2024). https://doi.org/10.1038/s41467-024-53260-y