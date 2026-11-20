options(stringAsFactors = F)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggupset)

c0 = function(x)
    cat(paste(x,"\n"))
t0 = function(x)
    cat("\t",paste(x,"\n"))
p0 = paste0

ds_dir = "diffsplice_bin/"
tissues = c("Cortex_Frontal","Cortex_Motor","Spinal_Cord_Cervical","Spinal_Cord_Lumbar","Cerebellum")
annots = read.table("all_tested_annots_ge30_exmatch.tsv",header=T,sep="\t",quote="",row.names=1)


get_sig = function(tissue) {
    sig = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_cluster_significance.txt"),
                     header=TRUE,sep="\t")
    sig = sig[which(sig$status=="Success"),]
    sig = sig[order(sig$p),]
    rownames(sig) = str_replace(sig$cluster,"chr.*:","")
    t0(p0("Tested clusters: ",nrow(sig)))
    sig = sig[which(sig$p.adjust<0.05),]
    t0(p0("Significant clusters: ",nrow(sig)))
    return(sig)
}
get_eff = function(tissue,sig) {
    eff = read.table(paste0(ds_dir,tissue,"/leafcutter_ds_effect_sizes.txt"),
                 header=TRUE,sep="\t",check.names=FALSE)
    eff$clu = str_replace(eff$intron,".*clu","clu")
    eff = eff[which(eff$clu %in% rownames(sig)),]
    eff = eff[order(abs(eff$deltapsi),decreasing=TRUE),]
    eff$clu_p.adjust = sig[eff$clu,]$p.adjust
    return(eff)
}
append_max = function(sig,eff) {
    sig$max_dpsi = -1
    for (i in 1:nrow(sig)) {
        sig$max_dpsi[i] = max(abs(eff[which(eff$clu==rownames(sig)[i]),]$deltapsi))
    }
    return(sig)
}

all_res = list()
all_eff = list()
sig_count = list()
for (tissue in tissues) {
    c0(tissue)
    groups = read.table(paste0(ds_dir,tissue,"/lc_groups.txt"),header=F,sep="\t")
    t0(p0("ALS n=",table(groups$V2)[1]))
    t0(p0("Control n=",table(groups$V2)[2]))
    sig = get_sig(tissue)
    eff = get_eff(tissue,sig)
    t0(p0("Significant clusters (∆PSI>0.01): ",length(which(rownames(sig) %in% unique(eff[which(abs(eff$deltapsi)>0.01),]$clu)))))
    #t0(p0("Significant clusters: ",nrow(sig)))
    t0(p0("Top genes: ",p0(head(sig$genes),collapse=";")))
    sig = append_max(sig,eff)
    sig$tissue = tissue
    all_res[[tissue]] = sig
    eff$tissue = tissue
    all_eff[[tissue]] = eff
    sig_count[[tissue]] = cbind.data.frame(tissue=tissue,sample_n=nrow(groups),sig_clus=nrow(sig))
    c0("")
}
all_res = do.call(rbind.data.frame,all_res)
sig_count = do.call(rbind.data.frame,sig_count)
all_eff = do.call(rbind.data.frame,all_eff)

all_eff$tsg_gene_name = annots[match(all_eff$intron,rownames(annots)),]$gene_name
all_eff$tsg_junc_match_rank = annots[match(all_eff$intron,rownames(annots)),]$junc_match_rank
all_eff$tsg_gene_id = annots[match(all_eff$intron,rownames(annots)),]$gene_id
all_eff$tsg_overlaps_cds = annots[match(all_eff$intron,rownames(annots)),]$overlaps_cds=="True"
all_res$cluster = str_replace_all(all_res$cluster,".*:","")
all_res = all_res[,colnames(all_res)[which(!colnames(all_res) %in% c("status","clu","genes"))]]

all_eff$deltapsi = all_eff$ALS-all_eff$Control

clus = cbind.data.frame(clu=unique(all_res$cluster),genes="",ov_cds=F)
for (i in 1:nrow(clus)) {
    clus$gene_names[i] = paste0(unique(unlist(str_split(all_eff[which(all_eff$clu==clus$clu[i]),]$tsg_gene_name,","))),collapse=",")
    clus$gene_ids[i] = paste0(unique(unlist(str_split(all_eff[which(all_eff$clu==clus$clu[i]),]$tsg_gene_id,","))),collapse=",")
    clus$ov_cds[i] = any(all_eff[which(all_eff$clu==clus$clu[i]),]$tsg_overlaps_cds)
}
rownames(clus) = clus$clu
all_res$tsg_gene_name = clus[all_res$cluster,]$gene_names
all_res$tsg_gene_id = clus[all_res$cluster,]$gene_ids
all_res$tsg_overlaps_cds = clus[all_res$cluster,]$ov_cds

write.table(all_res,"all_res.tsv",col.names=T,row.names=F,sep="\t")
write.table(all_res,"all_eff.tsv",col.names=T,row.names=F,sep="\t")