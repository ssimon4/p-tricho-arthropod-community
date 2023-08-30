rm(list=ls())
cat("\014")

#setwd("~/Users/harichhetri/Desktop/..")
library(lme4)

#name files:
tfam<-read.table("gatk_882_WG_genotypes_biallelic_snps_VQSR_tranches99_maf005_012.tfam")
Bnames<-read.table('id_to_name_882.txt');colnames(Bnames)<-c("JG1","Geno","JGIname")

#Phenotype files:
Ins<-read.table('clats_richness.txt')
#reformatting for only the traits of interest, renaming columns to match
Inspheno<-Ins[,c(1,2)];colnames(Inspheno)<-c("Geno","clats_richness")
#Creating merged phenotype matrix:
all.pheno<-merge(Inspheno,Bnames,all.y=T)
attach(all.pheno)

#creating the phenotype file in the correct order, same as the tfam file:
PMAT<-matrix(NA,1,1)
for(j in 1:nrow(tfam)){
	Pn<-which(all.pheno$JGIname==as.character(tfam[j,1]))
	Pmat<-all.pheno[Pn,2]
	PMAT<-rbind(PMAT,Pmat)
}
PFMAT_3col<-cbind(as.character(tfam[,1]),as.character(tfam[,1]),PMAT[-1,1])
PFMAT<- PFMAT_3col[,c(-1,-2)]
write.table(PFMAT,paste(colnames(Inspheno)[2],"_Pheno.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE)

