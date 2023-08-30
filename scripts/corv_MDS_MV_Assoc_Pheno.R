rm(list=ls())
cat("\014")

#setwd("~/Users/harichhetri/Desktop/..")
library(lme4)

#name files:
tfam<-read.table("gatk_882_WG_genotypes_biallelic_snps_VQSR_tranches99_maf005_012.tfam")
Bnames<-read.table('id_to_name_882.txt');colnames(Bnames)<-c("JG1","Geno","JGIname")

#Phenotype files:
Pheno<-read.table('corv_MDS_MV_outlier_removed_MADmethod.txt')
colnames(Pheno)<-c("Geno","MDS1", "MDS2", "MDS3", "MDS4", "MDS5")
#Creating merged phenotype and PC matrix:
PhPC<-merge(Pheno,Bnames,all.y=T)
attach(PhPC)

all_pheno<- PhPC[c(-1, -8)]
all_pheno<- all_pheno[c(6,1:5)]
#Multivariate linear regression for gemma covariates

all_pheno<- all_pheno[order(JGIname),]
assoc.all_pheno<- all_pheno[2:6]

write.table(assoc.all_pheno, "corv_MDS_MV_Pheno.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

