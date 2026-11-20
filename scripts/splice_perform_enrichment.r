options(stringAsFactors = F)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggupset)

library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

annots = read.table("all_tested_annots_ge30_exmatch.tsv",header=T,row.names=1,sep="\t",quote = "")
annots = annots[which(!grepl(",",annots$gene_id)),]

tested_de = unique(read.table("../als_diffsplice/deg/Cerebellum_als_vs_control.csv",header=T,sep=",")[,'X'])
tested_ds = unique(unlist(str_split(annots$gene_id,",")))
tested_ds = tested_ds[which(tested_ds %in% tested_de)]

to_entrez = function(x) { 
    return(unique(bitr(unique(str_replace_all(x,"\\..*","")),fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID))
}
trim_ensg = function(genes) {
    return(str_replace_all(genes,"\\..*",""))
}


ds_dir = "diffsplice_bin/"
tissues = c("Cortex_Frontal","Cortex_Motor","Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Cerebellum")

GO_TYPE="BP"
all_out = list()
for (tissue in tissues) {
    message(tissue)
    sig = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_cluster_significance.txt"),header=T,sep="\t")
    sig = sig[which(sig$status=="Success"),]
    sig$clu = str_replace_all(sig$clu,".*:","")
    eff = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_effect_sizes.txt"),header=T,sep="\t")
    eff$clu = str_replace_all(eff$intron,".*:","")
    tested_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    tested_ds_genes_tissue = tested_ds_genes_tissue[which(tested_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    sig = sig[which(sig$p.adjust<0.05),]
    eff = eff[which(eff$clu %in% sig$clu),]

    sig_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    sig_ds_genes_tissue = sig_ds_genes_tissue[which(sig_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    kk <- enrichGO(gene = trim_ensg(sig_ds_genes_tissue), 
                    universe = trim_ensg(tested_ds_genes_tissue),
                    OrgDb = org.Hs.eg.db,
                     keyType="ENSEMBL",
                     ont=GO_TYPE,
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "fdr",
                    readable=TRUE)
    if (nrow(kk)>0) {
        kk = cbind.data.frame(contrast="als_v_control",tissue=tissue,kk)
        all_out[[tissue]] = kk
    }
}
splice_enrich_go_bp = do.call(rbind.data.frame, all_out)


all_out = list()
for (tissue in tissues) {
    message(tissue)
    sig = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_cluster_significance.txt"),header=T,sep="\t")
    sig = sig[which(sig$status=="Success"),]
    sig$clu = str_replace_all(sig$clu,".*:","")
    eff = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_effect_sizes.txt"),header=T,sep="\t")
    eff$clu = str_replace_all(eff$intron,".*:","")
    tested_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    tested_ds_genes_tissue = tested_ds_genes_tissue[which(tested_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    sig = sig[which(sig$p.adjust<0.05),]
    eff = eff[which(eff$clu %in% sig$clu),]

    sig_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    sig_ds_genes_tissue = sig_ds_genes_tissue[which(sig_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    kk <- enrichKEGG(
                      gene          = to_entrez(sig_ds_genes_tissue),   # Significantly differentially expressed genes (as Entrez IDs)
                      universe      = to_entrez(tested_ds_genes_tissue), # All tested genes (as Entrez IDs)
                      organism      = "hsa",                            # Organism, e.g., "hsa" for human, "mmu" for mouse
                      pvalueCutoff  = 0.05,                             # P-value threshold for significance
                      pAdjustMethod = "fdr"                           # Multiple testing correction method (e.g., "fdr" for False Discovery Rate)
                    )
    if (nrow(kk)>0) {
        kk = cbind.data.frame(contrast="als_v_control",tissue=tissue,kk)
        all_out[[tissue]] = kk
    }
}
splice_enrich_kegg = do.call(rbind.data.frame, all_out)


all_out = list()
for (tissue in tissues) {
    message(tissue)
    sig = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_cluster_significance.txt"),header=T,sep="\t")
    sig = sig[which(sig$status=="Success"),]
    sig$clu = str_replace_all(sig$clu,".*:","")
    eff = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_effect_sizes.txt"),header=T,sep="\t")
    eff$clu = str_replace_all(eff$intron,".*:","")
    tested_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    tested_ds_genes_tissue = tested_ds_genes_tissue[which(tested_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    sig = sig[which(sig$p.adjust<0.05),]
    eff = eff[which(eff$clu %in% sig$clu),]

    sig_ds_genes_tissue = unique(unlist(str_split(annots[eff$intron,]$gene_id,",")))
    sig_ds_genes_tissue = sig_ds_genes_tissue[which(sig_ds_genes_tissue %in% tested_de)] # LIMIT TO DE

    kk <- enrichPathway(gene = to_entrez(sig_ds_genes_tissue), 
                    universe = to_entrez(tested_ds_genes_tissue),
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "fdr",
                    readable=TRUE)
    if (nrow(kk)>0) {
        kk = cbind.data.frame(contrast="als_v_control",tissue=tissue,kk)
        all_out[[tissue]] = kk
    }
}
splice_enrich_reactome = do.call(rbind.data.frame, all_out)

write.table(splice_enrich_reactome,"splice_enrich_reactome.tsv",col.names=T,row.names=F,sep="\t")
write.table(splice_enrich_kegg,"splice_enrich_kegg.tsv",col.names=T,row.names=F,sep="\t")
write.table(splice_enrich_go_bp,"splice_enrich_go_bp.tsv",col.names=T,row.names=F,sep="\t")
