# WGS Variant Calling Worfklows
This document outlines the commands used for processing the whole-genome data at the New York Genome Center.

---

## Overview

1. SNV & INDEL calling using GATK Best Practices.
2. Functional annotations using VEP.
3. Short Tandem Repeat (STR) genotyping using ExpansionHunter.
4. Structural Variant (SV) calling & genotyping:  
    4.1. Manta<br>
    4.2. Melt<br>
    4.3. Absinthe<br> 
    4.4. Paragraph<br>
    4.5. PanGenie<br>

---

## 1. SNV & INDEL calling using GATK Best Practices

### Alignment, post-processing, and SNV and INDEL calling
Alignment and post-processing are performed as outlined by the Center for Common Disease Genomics project: [https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md](https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md) .

### Programs and reference data
The data was aligned to the reference genome using the following programs and reference datasets:

1. [BWA-MEM](https://github.com/lh3/bwa/blob/master/bwakit/README.md)
2. [Samtools-1.3.1](http://www.htslib.org/doc/samtools-1.3.1.html)
3. [Picard-2.4.1](https://github.com/broadinstitute/picard/releases/tag/2.4.1)
3. [GATK-3.5-0](https://github.com/broadgsa/gatk-protected/tree/3.5)
4. Resource files can be obtained here: [gs://gcp-public-data--broad-references/hg38/v0](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0) .

### Reference genome: GRCh38 with alternative sequences, plus decoys and HLA
The reference genome that the data was aligned to can be obtained here: [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)

### Command lines

1. Alignment at lane level

```bash
bwa mem \
  -Y \
  -K 100000000 \
  -t 16 \
  -R $rg_string \
  $reference_fasta \
  $fastq_file(1) \
  $fastq_file(2) | samtools view -Shb -o $bam_aligned -
```

2. Fix mate information in the BAM

```bash       
java $jvm_args -jar picard.jar \
  FixMateInformation \
  MAX_RECORDS_IN_RAM=2000000 \
  VALIDATION_STRINGENCY=SILENT \
  ADD_MATE_CIGAR=True \
  ASSUME_SORTED=true \
  I=$bam_aligned \
  O=$bam_fixmate
```

3. Merging lane-level BAM files to generate a single sample-level BAM file

```bash
java $jvm_args -jar picard.jar \
  MergeSamFiles \
  USE_THREADING=true \
  MAX_RECORDS_IN_RAM=2000000 \
  VALIDATION_STRINGENCY=SILENT \
  SORT_ORDER=queryname \
  INPUT=$bam_fixmate_1 \
  INPUT=$bam_fixmate_2 \
  OUTPUT=$bam_merged
```

4. Mark duplicates and coordinate-sort the BAM

```bash
java $jvm_args -jar picard.jar \
  MarkDuplicates \
  MAX_RECORDS_IN_RAM=2000000 \
  VALIDATION_STRINGENCY=SILENT \
  M=$dedup_metrics \
  I=$bam_merged \
  O=$bam_dedup   
      
```
```bash
java $jvm_args -jar picard.jar \
  SortSam \
  MAX_RECORDS_IN_RAM=2000000 \
  VALIDATION_STRINGENCY=SILENT \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true \
  I=$bam_dedup \
  O=$bam_sorted
```

5. Recalibrate base quality scores using known SNPs and generate a BAM with recalibrated base qualities

```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -downsample_to_fraction 0.1 \
  -nct 4 \
  --preserve_qscores_less_than 6 \
  -L $autosomes \
  -R $reference_fasta \
  -o $recal_data.table \
  -I $bam_sorted \
  -knownSites $known_snps_from_dbSNP138 \
  -knownSites $known_indels \
  -knownSites $known_indels_from_mills_1000genomes
```
```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T PrintReads \
  -nct 4 \
  --disable_indel_quals \
  --preserve_qscores_less_than 6 \
  -SQQ 10 \
  -SQQ 20 \
  -SQQ 30 \
  -rf BadCigar \
  -R $reference_fasta \
  -o $bam_recalibrated \
  -I $bam_sorted \
  -BQSR $recal_data.table
```

6. Creating CRAM file

```bash
samtools view \
  -C \
  -T $reference_fasta \
  -o $cram \
  $bam_recalibrated
```
```bash
samtools index $cram
```

7. SNV/INDEL discovery using HaplotypeCaller in genomic-VCF (GVCF) mode

For variant discovery, we used HaplotypeCaller in GVCF mode (Poplin et al., 2017) with sex-dependent ploidy settings on chromosome X and Y. Specifically, variant discovery on chrX was performed using diploid settings in females, diploid settings on PAR regions in males, and haploid settings on non-PAR regions in males. Variant discovery on chrY was performed with haploid settings in males and skipped in females.

```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  --genotyping_mode DISCOVERY \
  -A AlleleBalanceBySample \
  -A DepthPerAlleleBySample \
  -A DepthPerSampleHC \
  -A InbreedingCoeff \
  -A MappingQualityZeroBySample \
  -A StrandBiasBySample \
  -A Coverage \
  -A FisherStrand \
  -A HaplotypeScore \
  -A MappingQualityRankSumTest \
  -A MappingQualityZero \
  -A QualByDepth \
  -A RMSMappingQuality \
  -A ReadPosRankSumTest \
  -A VariantType \
  -l INFO \
  --emitRefConfidence GVCF \
  -rf BadCigar \
  --variant_index_parameter 128000 \
  --variant_index_type LINEAR \
  -R $reference_fasta \
  -nct 1 \
  -L $interval \
  -ploidy $ploidy \
  -I $bam_recalibrated \
  -o $gvcf
```
8. Combining GVCF files (executed in batches of 200 samples)

```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T CombineGVCFs \
  -R $reference_fasta \
  -L $interval \
  --variant $gvcf_list \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -o $combined_gvcf
```

9. Joint Genotyping across all samples

```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T GenotypeGVCFs \
  -R $reference_fasta \
  -nt 5 \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L $interval \
  --variant $combined_gvcf_list \
  -o $vcf_genotyped \
```

10. Variant Quality Score Recalibration (VQSR) to assign FILTER status

```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R $reference_fasta \
  -nt 5 \
  -input $vcf_genotyped \
  -mode SNP \
  -recalFile $vqsr_snp.recal \
  -tranchesFile $vqsr_snp.tranches \
  -rscriptFile $vqsr_snp_plots.R \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 $kg_omni \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kg_snps \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD \
  -an MQ \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -an SOR \
  -an DP \
  -tranche 100.0 \
  -tranche 99.8 \
  -tranche 99.6 \
  -tranche 99.4 \
  -tranche 99.2 \
  -tranche 99.0 \
  -tranche 95.0 \
  -tranche 90.0 \
```
```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R $reference_fasta \
  -nt 5 \
  -input $vcf_genotyped \
  -mode INDEL \
  -recalFile $recalibrate_indel.recal \
  -tranchesFile $recalibrate_indel.tranches \
  -rscriptFile $recalibrate_indel_plots.R \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 $kg_mills \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
  -an QD \
  -an FS \
  -an ReadPosRankSum \
  -an MQRankSum \
  -an SOR \
  -an DP \
  -tranche 100.0 \
  -tranche 99.0 \
  -tranche 95.0 \
  -tranche 92.0 \
  -tranche 90.0 \
  --maxGaussians 4
```
```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -R $reference_fasta \
  -nt 5 \
  -input $vcf_genotyped \
  -mode SNP \
  --ts_filter_level 99.80 \
  -recalFile $recalibrate_SNP.recal \
  -tranchesFile $recalibrate_SNP.tranches \
  -o $vcf_recalibrated_snp \
```
```bash
java $jvm_args -jar GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -R $reference_fasta \
  -nt 5 \
  -input $vcf_recalibrated_snp \
  -mode INDEL \
  --ts_filter_level 99.0 \
  -recalFile $recalibrate_INDEL.recal \
  -tranchesFile $recalibrate_INDEL.tranches \
  -o $vcf_recalibrated_snp_indel \
```

## 2. Functional annotations using Variant Effect Predictor (VEP)

- VEP v114.0 was run using Singularity v4.2.2.
- fastaReference="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" (provided with VEP installation).
- clinvarVCF: clinvar resource file was downloaded from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/ (fileDate=2025-05-04)
- num_forks=4

```bash
singularity exec vep.sif \
  vep --dir $HOME/vep_data \
  --cache \
  --offline \
  --exclude_predicted \
  --fasta $fastaReference \
  --gencode_basic \
  --sift p \
  --polyphen p \
  --hgvs \
  --symbol \
  --numbers \
  --domains \
  --regulatory \
  --tsl \
  --gene_phenotype \
  --pubmed \
  --nearest symbol \
  --pick_allele_gene \
  --af_gnomade \
  --af_gnomadg \
  --vcf \
  --compress_output bgzip \
  --force_overwrite \
  --fork $num_forks \
  --input_file $vcf \
  --output_file $prefix.vep_v114.vcf.gz \
  --custom file=$clinvarVCF,format=vcf,short_name=clinvar_20250504,type=exact,overlap_cutoff=0,fields=DBVARID%ALLELEID%CLNDN%CLNDISDB%MC%CLNSIG%CLNSIGCONF%CLNREVSTAT%ORIGIN
```

## 3. Short Tandem Repeat (STR) genotyping using ExpansionHunter

### Programs and reference data
- [ExpansionHunter 5.0](https://github.com/Illumina/ExpansionHunter/releases/tag/v5.0.0)
- [VariantCatalog](https://github.com/nygenome/nygc-germline-pipeline-readmes/blob/master/EHv5_variant_catalog.json)
- [Reference Genome](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
   
```bash
$EH \
  --reads $BAM \
  --reference $REF \
  --variant-catalog $JSON \
  --output-prefix $S.eh \
  --sex $sex
```

## 4. Structural Variant (SV) calling & genotyping

### 4.1 Manta
- [Manta 1.5.0](https://github.com/Illumina/manta/releases/tag/v1.5.0)
- [Reference Genome](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
- Threads = 8
- Mem = 16

```bash
/manta-1.5.0/bin/configManta.py \
  --bam $bam \
  --referenceFasta $ref \
  --runDir $sample

execute=$sample/runWorkflow.py

$execute \
  --mode local \
  --jobs $threads \
  --memGb $mem
```

### 4.2. Melt

- JAVAVER = jdk-1.8.0.45
- JAR = MELTv2.2.2
- [Reference Genome](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)
- meiType = ALU, LINE1, HERVK, SVA
- zipPath = /MELTv2.2.2/me_refs/Hg38/${meiType}_MELT.zip
- BOWTIE = bowtie2-2.5.1
- GENES = MELTv2.2.2/add_bed_files/Hg38/Hg38.genes.bed

Preprocessing

```bash
$JAVAVER -Xmx20G -jar $JAR \
  Preprocess \
  -bamfile $cram \
  -h $REF;
```
MEI calling 

```bash
$JAVAVER -Xmx20G -jar $JAR \
  IndivAnalysis \
  -r 150 \
  -bamfile $cram \
  -w $outdir/temp/$meiType \
  -t $zipPath \
  -h $REF \
  -bowtie $BOWTIE;
```

Joint calling

```bash
$JAVAVER -Xmx100G -jar $JAR \
  GroupAnalysis \
  -r 150 \
  -discoverydir $outdir/temp/$meiType \
  -w $outdir/temp/$meiType \
  -n $GENES \
  -t $zipPath \
  -h $REF;
```

Genotyping

```bash
$JAVAVER -Xmx100G -jar $JAR \
  Genotype \
  -bamfile $cram \
  -p $outdir/temp/$meitype \
  -w $outdir/temp/$meitype \
  -t $ZIP \
  -h $REF;
```

VCF output

```bash
$JAVAVER -Xmx20G -jar $JAR \
  MakeVCF \
  -ac \
  -genotypingdir $outdir/temp/$meiType \
  -h $REF \
  -t $zipPath \
  -w $outdir/temp/$meiType \
  -p $outdir/temp/$meiType;
```

### 4.3. Absinthe

Resources - reference and alignment resource files

```bash
  absinthe \
  -s $sample \
  -c $cram \
  -r $ref \
  -R /resources \
  -o /out \
  -d $depth \
  -l $readlen \
  -g $sex \
  -t $THREADS

```

### 4.4. Paragraph
MAXDEPTH = 5 * sample coverage

Conversion of VCF to JSON
```bash
convertVCF2JSON.py -i $vcf \
  -r ${REF} \
  -o /out \
  --threads ${THREADS}
```
Genotyping
```bash
multigrmpy_tmpjson.py -i $vcfJson \
  -m $manifest \
  -r ${REF} \
  -o $out \
  --threads ${THREADS} \
  -M ${MAXDEPTH} \
  --scratch-dir $tmp
```

### 4.5. PanGenie
- VCF - [HGSVC](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/PanGenie_PAV-panel/20210311_pav-panel-freeze4.vcf.gz). 
- biallelicVCF - [HGSVC](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200825_HHU_genotyping-vcf-freeze3/20200825_pav-all-samples-freeze3-final.vcf.gz). 

Genotyping
```bash
  PanGenie \
  -i $fq \
  -r $ref \
  -v $vcf \
  -o /out \
  -s $sample \
  -t $THREADS \
  -j $THREADS
```
Conversion to biallelic
```bash
zcat /out/${sample}_genotyping.vcf.gz \
  | convert-to-biallelic.py \
    $biallelicVCF
```
