# ---
# title: STR QTL analysis for clinical STRs
# author: Yebin Kim
# date: 2025-08-12
# ---

import pandas as pd
import numpy as np
import os
import re
import tensorqtl
from tensorqtl import cis


## load variant df
eh_table_als = pd.read_csv("cgnd_n4746_EH_STRMerge_WO_NonConsentedSamps_WOQcOutliers_removedXY_filtered.tsv",sep='\t')
eh_table_als = eh_table_als.dropna(axis=1)
als_header = pd.read_csv("header_4746.txt",sep='\t',header=None,index_col=None)
als_header = list(als_header.loc[0,:])
eh_table_als.columns = als_header
eh_table_als = eh_table_als[eh_table_als.FILTER == 'PASS'] # passed dumpSTR QC only

variant_df = eh_table_als[['#CHROM','POS']]
variant_df.index = eh_table_als['REPID']
variant_df['#CHROM'] = variant_df['#CHROM'].str.replace("chr","")
variant_df.columns = ['chrom','pos']
variant_df.chrom =  variant_df.chrom.astype(int)
variant_df.pos =  variant_df.pos.astype(int)
variant_df = variant_df.sort_values(['chrom','pos'])
variant_df.chrom = variant_df.chrom.astype(str)


# load genotype df
max_allele_table = pd.read_csv("max_allele_table.tsv",sep='\t',header=0,index_col=0) 
max_allele_table = max_allele_table.apply(pd.to_numeric,errors='coerce') 

min_allele_table = pd.read_csv("min_allele_table.tsv",sep='\t',header=0,index_col=0) 
min_allele_table = min_allele_table.apply(pd.to_numeric,errors='coerce')

# include pcr free samples only
pcrfree_samples = pd.read_csv("wgs_pcr_free_samples.csv",index_col=None,header=None)
pcrfree_samples = list(pcrfree_samples[0]) 


def run_QTL(pcrfree_samples, max_allele_table, min_allele_table, variant_df, outdir, group):

    g_bed_file =f"geneexp_pheno/NYGC_{group}_pheno.regressed.harmonised.bed" # gene expression matrix
    s_bed_file =f"splicing_pheno/NYGC_{group.lower()}_pheno.regressed.bed" # splicing matrix

    pheno_samples = pd.read_csv(g_bed_file, header=0, index_col=None, nrows=1,sep='\t')
    overlapping_samples = list(set(pheno_samples.columns) & set(pcrfree_samples))

    # ----------------- genotype data ---------------------
    cate_max_allele_table = max_allele_table[overlapping_samples].copy()
    cate_min_allele_table = min_allele_table[overlapping_samples].copy()
    cate_max_allele_table = cate_max_allele_table[np.var(cate_max_allele_table,axis=1)>=1]
    cate_min_allele_table = cate_min_allele_table[cate_min_allele_table.index.isin(list(cate_max_allele_table.index))]
    
    variant_df = variant_df.loc[list(cate_max_allele_table.index)] # match the variants

    threshold_1 = []  
    threshold_2 = []     
    perc_0 = []  
    perc_1 = []  
    perc_2 = []    
    
    for i in range(cate_max_allele_table.shape[0]): # for each variant
        variant = cate_max_allele_table.index[i]
        mean = np.nanmean(cate_min_allele_table.loc[variant])
        window = np.std(cate_max_allele_table.loc[variant]) 
 
        left_bound = mean
        right_bound = mean+window

        threshold_1 += [left_bound]
        threshold_2 += [right_bound]

        count_0 = 0
        count_1 = 0
        count_2 = 0    

        for j in range(cate_max_allele_table.shape[1]): # for each sample
            value = cate_max_allele_table.iat[i, j]
            if pd.isna(value): 
                continue  

            if value >= right_bound:
                cate_max_allele_table.iat[i, j] = 2
                count_2 += 1

            elif value <= left_bound:
                cate_max_allele_table.iat[i, j] = 0
                count_0 += 1

            elif left_bound < value < right_bound:
                cate_max_allele_table.iat[i, j] = 1
                count_1 += 1

        perc_0 += [round(count_0/max_allele_table.shape[1]*100,1)]
        perc_1 += [round(count_1/max_allele_table.shape[1]*100,1)]
        perc_2 += [round(count_2/max_allele_table.shape[1]*100,1)]

    cate_max_allele_table = cate_max_allele_table.fillna(-9) # for tensorQTL imputation

    # ----------------- covariate data ------------------------
    covariate_df = pd.read_csv(f"pcrfree_pca_cov_{group}.tsv",sep='\t',header=0,index_col=0)
    covariate_df = covariate_df.loc[overlapping_samples]

    # ----------------- phenotype data ------------------------
    # (1) eQTL
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(g_bed_file)
    phenotype_pos_df.start = phenotype_pos_df.start.astype(int)
    phenotype_pos_df.end = phenotype_pos_df.end.astype(int)
    phenotype_pos_df.chr = phenotype_pos_df.chr.str.replace("chr","")
    phenotype_df = phenotype_df[overlapping_samples]

    egenes_df = cis.map_cis(cate_max_allele_table, variant_df, phenotype_df ,phenotype_pos_df, covariate_df ,nperm=10000)
    egenes_df.to_csv(f"{outdir}/eGenes_eQTL_clinicalSTR_{group}.tsv", sep='\t',header=True, index=True)


    # (2) sQTL
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(s_bed_file)
    phenotype_pos_df.start = phenotype_pos_df.start.astype(int)
    phenotype_pos_df.end = phenotype_pos_df.end.astype(int)
    phenotype_pos_df.chr = phenotype_pos_df.chr.str.replace("chr","")
    phenotype_df = phenotype_df[overlapping_samples]

    egenes_df = cis.map_cis(cate_max_allele_table, variant_df, phenotype_df ,phenotype_pos_df, covariate_df ,nperm=10000)
    egenes_df.to_csv(f"{outdir}/eGenes_sQTL_clinicalSTR_{group}.tsv", sep='\t',header=True, index=True)


directory = "/path/to/data"
output_directory = "/path/to/output"
pattern = re.compile(r'NYGC_(.*?)_pheno\.regressed\.harmonised\.tsv')

for filename in os.listdir(directory):
    if filename.endswith(".tsv"):
        match = pattern.search(filename)
        if match:
            group_name = match.group(1)
            file_path = os.path.join(directory, filename)
            run_QTL(pcrfree_samples, max_allele_table, min_allele_table, variant_df, output_directory, group_name)
