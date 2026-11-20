options(stringAsFactors = F)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggupset)

covs = read.table("metadata_updated_sep4_2024_with_merged_mc.csv",
                 header=T,sep=",")
covs = covs[order(covs$RIN,decreasing = T),]
covs = covs[order(covs$Platform,decreasing = T),]
covs$Site.Specimen.Collected = str_replace_all(covs$Site.Specimen.Collected," ","_")
covs$Prep = str_replace_all(covs$Prep," ","_")


TISSUE = "Cerebellum"
tmp = covs[which(covs$Sample.Source=="Cerebellum"),]
dim(tmp)
tmp = tmp[which(!duplicated(tmp$ExternalSubjectId)),]
dim(tmp)
tmp = tmp[,c("RNA","Subject.Group","Prep","MEAN_READ_LENGTH","Site.Specimen.Collected")]
write.table(tmp,"diffsplice_bin/1_Cerebellum/lc_groups.txt",col.names=F,row.names=F,sep=" ",quote=F)

TISSUE = "Spinal_Cord_Lumbar"
tmp = covs[which(covs$Sample.Source==TISSUE),]
dim(tmp)
tmp = tmp[which(!duplicated(tmp$ExternalSubjectId)),]
dim(tmp)
tmp = tmp[,c("RNA","Subject.Group","Site.Specimen.Collected","Prep","PCT_CHIMERAS")]
head(tmp,3)
table(is.na(tmp))
write.table(tmp,paste0("diffsplice_bin/2_",TISSUE,"/lc_groups.txt"),col.names=F,row.names=F,sep=" ",quote=F)

TISSUE = "Spinal_Cord_Cervical"
tmp = covs[which(covs$Sample.Source==TISSUE),]
dim(tmp)
tmp = tmp[which(!duplicated(tmp$ExternalSubjectId)),]
dim(tmp)
tmp = tmp[,c("RNA","Subject.Group","Prep","MEAN_READ_LENGTH","Site.Specimen.Collected")]
head(tmp,3)
table(is.na(tmp))
write.table(tmp,paste0("diffsplice_bin/3_",TISSUE,"/lc_groups.txt"),col.names=F,row.names=F,sep=" ",quote=F)

TISSUE = "Cortex_Frontal"
tmp = covs[which(covs$Sample.Source==TISSUE),]
dim(tmp)
tmp = tmp[which(!duplicated(tmp$ExternalSubjectId)),]
dim(tmp)
tmp = tmp[,c("RNA","Subject.Group","Prep","MEAN_READ_LENGTH","PCT_CHIMERAS")]
head(tmp,3)
table(is.na(tmp))
write.table(tmp,paste0("diffsplice_bin/5_",TISSUE,"/lc_groups.txt"),col.names=F,row.names=F,sep=" ",quote=F)

TISSUE = "Cortex_Motor"
tmp = covs[which(covs$Sample.Source==TISSUE),]
dim(tmp)
tmp = tmp[which(!duplicated(tmp$ExternalSubjectId)),]
dim(tmp)
tmp = tmp[,c("RNA","Subject.Group","Prep","MEAN_READ_LENGTH","PCT_CHIMERAS")]
head(tmp,3)
table(is.na(tmp))
write.table(tmp,paste0("diffsplice_bin/7_",TISSUE,"/lc_groups.txt"),col.names=F,row.names=F,sep=" ",quote=F)




