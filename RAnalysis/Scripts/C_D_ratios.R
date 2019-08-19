#Title: PGA QPCR Symbiont ID Data
#Author: HM Putnam
#Date Last Modified: 20170504
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
##### Include Versions of libraries


#Required Data files

# Set Working Directory:
setwd("~/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/") #set working
mainDir<-"~/MyProjects/Mcap_PGA_TGA/RAnalysis/" #set main directory

# =================================================================================================
# DATA PREPARATION
# =================================================================================================
# • Load libraries --------------------------------------------------------------------------------
library(lme4); library(MASS); library(reshape2); library(lattice); library(lmerTest); library(lsmeans)
library(pbkrtest); library(scales); library(merTools); library(devtools); library(pBrackets)
## SPIDA package available at http://r-forge.r-project.org/projects/spida/
#system(paste("svn checkout svn://svn.r-forge.r-project.org/svnroot/spida/"))
#devtools::install("spida/pkg")
library(spida)

addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

# -------------------------------------------------------------------------------------------------
# • Load data -------------------------------------------------------------------------------------
# Use steponeR to import data and calculate S/H ratios
source("~/MyProjects/Mcap_PGA_TGA/RAnalysis/Scripts/QPCR.R")
#source(https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R) 
# Get list of plate files to read in

#file1 <- read.csv("~/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/qPCR/20170301_TGA2016_qPCR_1.csv", sep=",", header=F, skip=7)
#file2 <- read.csv("~/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/qPCR/20170313_TGA2016_qPCR_2.csv", sep=",", header=F, skip=21)
#Mcap.plates <- rbind(file1, file2)


Mcap.plates <- list.files(path="~/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/qPCR", pattern="csv$", full.names=T)


Mcap <- steponeR(files=Mcap.plates, target.ratios=c("C.Mcap", "D.Mcap"),
                 fluor.norm=list(C=2.26827, D=0, Mcap=0.84815),
                 copy.number=list(C=33, D=3, Mcap=1),
                 ploidy=list(C=1, D=1, Mcap=2), 
                 extract=list(C=0.813, D=0.813, Mcap=0.982))
Mcap <- Mcap$result
Mcap <- Mcap[-c(16,17),] #remove NTC and POS control
Mcap$D.CT.mean[is.na(Mcap$D.CT.mean)] <- 40 #set no amplification to max Ct
Mcap$C.CT.mean[is.na(Mcap$C.CT.mean)] <- 40 #set no amplification to max Ct
Mcap$ratio <- 2^(Mcap$D.CT.mean-Mcap$C.CT.mean) #Calculate the C:D ratio
Mcap$prop.C <- Mcap$ratio/(Mcap$ratio+1) #calculate proportion of C
Mcap$prop.D <- 1-Mcap$prop.C #calculate proportion of D
Mcap$prop.C <-round(Mcap$prop.C, digits = 2) #round values to 2 decimal places
Mcap$prop.D <-round(Mcap$prop.D, digits = 2) #round values to 2 decimal places

props <- data.frame(Mcap$Sample.Name, Mcap$prop.C, Mcap$prop.D) #create a dataframe of sampleID and proportions
colnames(props) <- c("Sample.Name", "Prop.C", "Prop.D") #label columns
props$History <- c("Bleached", "Bleached", "Bleached", "Not Bleached", "Not Bleached", "Not Bleached", "Not Bleached", "Not Bleached", "Not Bleached", "Bleached", "Bleached", "Bleached","Bleached", "Bleached", "Bleached")
props$Treatment <-c("Initial","Final", "Final", "Initial","Final", "Final", "Initial","Final", "Final", "Initial","Final", "Final", "Initial","Final", "Final")
par(mar=c(10, 3, 1, 1))
boxplot(props$Prop.C ~props$Treatment*props$History, horizontal=F, las=2)

## Initial Clades

init <- subset(props, Treatment=="Initial")
boxplot(init$Prop.C ~init$History)

initial <- aggregate(Prop.C ~History, mean, data=props)
initial$Prop.D <- 1-initial$Prop.C
prop <- prop.table(initial[,2:3])
prop <- prop*2

barplot(as.matrix(t(prop)), xlab="Proportion C", col=c("black","gray"),
        legend = rownames(prop)) 

## Final Clades

fin <- subset(props, Treatment=="Final")
boxplot(fin$Prop.C ~fin$History)

final <- aggregate(Prop.C ~History, mean, data=fin)
final$Prop.D <- 1-final$Prop.C
propf <- prop.table(final[,2:3])
propf <- propf*2


jpeg("../Output/Init_Fin_Clades.jpg")
par(mfrow=c(1,2))
barplot(as.matrix(t(prop)), col=c("black","gray")) 
barplot(as.matrix(t(propf)), col=c("black","gray"),
        legend = colnames(propf)) 
dev.off()
