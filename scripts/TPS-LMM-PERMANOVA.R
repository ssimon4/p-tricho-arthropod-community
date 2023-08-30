#R-script for using TPS to adjust for micro-site variation, and calculating heritability.
library(lme4)
library(RLRsim)
library(spam)
library(fields)
library(msm)

setwd("C:/Users/ssimon4/Desktop/New Phytologist/RawData")

raw_data <- read.csv('ArthropodSurveyData.csv', header = T)

july_2012<-subset(raw_data, Year=='Jul-12')

Clatskanie<-subset(july_2012, Garden=='Clatskanie')

#TPS Corrections
#all_data<-read.csv('YourData.csv',header=T)
geno_trait<- subset(Clatskanie, select=c("Geno","Block", "Row", "Column", "Harmandia.Cecidomyid")) #Subsetting the required columns
attach(geno_trait)
Row<-as.numeric(geno_trait$Row)
Col<-as.numeric(geno_trait$Col)
#Rep<-as.factor(geno_trait$Rep)

jpeg(paste("TPS",colnames(geno_trait[5]),".jpg",sep=""), width=7.5, height=7.5, units="in", res=800)
par(mfrow=c(2,2))
plot(fitTps)
dev.off()

jpeg(paste("TPS_surface",colnames(geno_trait[5]),".jpg",sep=""), width=7.5, height=7.5, units="in", res=800)
out.p<-predictSurface(fitTps)
surface(out.p,type="C",xlab="Row",ylab="Col")
dev.off()

#fit Thin Plate Spline (TPS)
fitTps<-Tps(cbind(Row,Col),geno_trait[,5])
geno_trait<-cbind(geno_trait,resid(fitTps),resid(fitTps)+mean(geno_trait[,5]))
colnames(geno_trait)[6:7]<-c(paste('ResidTps_',colnames(geno_trait[5]),sep=""),paste('Tps_adjusted_',colnames(geno_trait[5]),sep=""))
#write.csv(geno_trait,'CorrectedValues.csv',row.names=FALSE)

clatskanie_2012_replicates<-read.csv('Clatskanie2012.csv', header = T)

#Linear mixed model 
H2mat<-matrix(NA,1,5);rownames(H2mat)<-colnames(clatskanie_2012_replicates[12]);colnames(H2mat)<-c("H2","GenoV", "ErrorV", "RLRatio","RLRTpvalue")
attach(clatskanie_2012_replicates)
#REML model for TPS corrected phenotypes 
fit<-lmer(clatskanie_2012_replicates[,12] ~ (1|Latitude) + (1|Geno))
fit2<-lmer(clatskanie_2012_replicates[,12] ~ Latitude + (1|Geno))

#Heritability calculation and significance test
#summary(fit)
Gvar<-VarCorr(fit)$Geno[1,1]
Evar<-sigma(fit)^2
H2<-Gvar/(Gvar+Evar)
H2mat[colnames(geno_trait[7]),1:3]<-c(H2,Gvar,Evar)  
LRTstats<-exactRLRT(fit2)
H2mat[colnames(geno_trait[7]),4]<-LRTstats$statistic
H2mat[colnames(geno_trait[7]),5]<-LRTstats$p.value

#R-script for running non-metric multidimensional scaling (NMDS) and PERMANOVA (ADONIS)
library(vegan);
library(labdsv);
library(ggplot2);
library(ggrepel);
library(plotly);
library(grid);
library(R.devices);
library(goeveg);
library(vegan3d);
library(processx);
library(plot3D);
library(car);
library(rgl)

clatskanie_2012_community<-read.table(file="Clatskanie2012Community.csv",header=TRUE,sep=",")
attach(clatskanie_2012_community)

NMDS5Clats<-metaMDS(clatskanie_2012_community[,11:27],k=5,distance="bray") # Dimensions = 5, Stress = 0.13

NMDSplot = data.frame(MDS1 = NMDS5Clats$points[,1], MDS2 = NMDS5Clats$points[,2], MDS3 = NMDS5Clats$points[,3], MDS4 = NMDS5Clats$points[,4], MDS5 = NMDS5Clats$points[,5])

write.table(file="Clatskanie2012NMDSPoints.csv", NMDSp, sep =",", row.names=TRUE)

shannon<-diversity(clatskanie_2012_community[,11:27],index="shannon")
write.table(file="ClatskanieShannon2012.csv", shannon, sep =",", row.names=TRUE)

jpg(file="Clats12-3dPlotFinal.jpg")
scatter3D(NMDSplot$MDS1, NMDSplot$MDS2,NMDSplot$MDS3, theta = 130, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed", colvar = NULL, col="darkred")
dev.off()


#PERMANOVA analysis of insect community in individual gardens
#Notes: genotypes that lack latitude information were removed from data prior to analysis

Clatskanie_Community_2012<-read.table(file="ClatskanieCommunity2012.csv.csv",header=TRUE,sep=",")

#PERMANOVA with more than one factor:
Clatskanie2012adonis<-adonis(clatskanie_2012_community[,11:27] ~ Latitude + Geno + as.factor(Block), data=clatskanie_2012_community)
Clatskanie2012adonis$aov.tab


#Run NMDS and ADONIS
all_gardens_2012<-read.table(file="AllGardens2012.csv",header=TRUE,sep=",")
attach(all_gardens_2012)

NMDS4All2012<-metaMDS(all_gardens_2012[,11:33],k=4,distance="bray") # Dimensions = 3, Stress = 0.10

NMDSAllPlot = data.frame(MDS1 = NMDS4All2012$points[,1], MDS2 = NMDS4All2012$points[,2], MDS3 = NMDS4All2012$points[,3])

pdf(file="AllGardens12-3dPlotFinal.pdf")
scatter3D(NMDSAllPlot$MDS1, NMDSAllPlot$MDS2,NMDSAllPlot$MDS3, theta = 130, bty = "g",
          pch = 20, cex = 2, ticktype = "detailed", colvar = NULL, col= all_gardens_2012$Color)
dev.off()

AllGardensAdonis<-adonis(all_gardens_2012[,11:33] ~ Latitude + Garden + Geno + as.factor(Block), data=all_gardens_2012)
AllGardensAdonis$aov.tab

#Generate interactive figure with plotly
fig <- plot_ly(NMDSAllPlot, x = NMDSAllPlot$MDS1, y = NMDSAllPlot$MDS2, z = NMDSAllPlot$MDS3, color = all_gardens_2012$Garden, colors = c('#2345d8', '#990d0d', '#ceae25'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'NMDS Axis 1'),
                                   yaxis = list(title = 'NMDS Axis 2'),
                                   zaxis = list(title = 'NMDS Axis 3')))

