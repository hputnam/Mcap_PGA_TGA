#Title: PGA Respirometry Data
#Author: HM Putnam
#Date Last Modified: 20160920
#See Readme file for details

rm(list=ls()) #clears workspace 

#Read in required libraries
##### Include Versions of libraries

library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")

#Required Data files

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/") #set working
mainDir<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/" #set main directory

#PHTOTSYNTHESIS
path.p<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/Respirometry/Photo" #the location of all your respirometry files
file.names<-list.files(path = path.p, pattern = "csv$") #list all csv file names in the folder
Photo.R <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Fragment.ID","Intercept", "µmol.L.sec")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data <-read.table(file.path(path.p,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data  <- Photo.Data [,c(9,16)] #subset columns of interest
  n<-dim(Photo.Data )[1] #identify length of data
  Photo.Data <-Photo.Data [60:(n-3),] #start at data point 1 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data )[1] #list length of trimmed data
  Photo.Data $sec <- 1:n #set seconds by one from start to finish of run
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/photo.thinning.pdf")
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  
  Photo.Data   <-  thinData(Photo.Data , by=10)$newData1 #thin data by every 10 points
  Photo.Data $sec <- as.numeric(rownames(Photo.Data )) #maintain numeric values for time
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  dev.off()
  
  Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $Value, alpha=0.2, 
                       method="pc", verbose=TRUE) 
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/photo.figs.pdf")
  plot(Regs)
  dev.off()
  
  Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[i,1] <- sub("_.*", "", file.names[i]) #stores the file name in the Date column
}

Photo.R

#RESPIRATION

path.r<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/Respirometry/Resp/" #the location of all your respirometry files
file.names.r<-list.files(path = path.r, pattern = "csv$")
Resp.R <- data.frame(matrix(NA, nrow=length(file.names.r), ncol=3, dimnames=list(file.names.r,c("Fragment.ID","Intercept", "µmol.L.sec")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names.r)) { # for every file in list start at the first and run this following function
  Resp.Data <-read.table(file.path(path.r,file.names.r[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Resp.Data  <- Resp.Data[,c(9,16)] #subset columns of interest
  n<-dim(Resp.Data)[1] #identify length of data
  Resp.Data <-Resp.Data[60:(n-3),] #start at data point 1 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Resp.Data)[1] #list length of trimmed data
  Resp.Data$sec <- 1:n #set seconds by one from start to finish of run
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/resp.thinning.pdf")
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Resp.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
  #set plot layout info
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Resp.Data$Value ~ Resp.Data$sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  
  Resp.Data  <-  thinData(Resp.Data, by=10)$newData1
  Resp.Data$sec <- as.numeric(rownames(Resp.Data))
  plot(Value ~ sec, data=Resp.Data, xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Resp.Data$Value ~ Resp.Data$sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  dev.off()
  
  Regs  <-  rankLocReg(xall=Resp.Data$sec, yall=Resp.Data$Value, alpha=0.2, 
                       method="pc", verbose=TRUE) 
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/resp.figs.pdf")
  plot(Regs)
  dev.off()
  
  Resp.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts them in the dataframe
  Resp.R[i,1] <- sub("_.*", "", file.names.r[i]) #stores the sample name in the Sample.ID column
}

Resp.R #view data

#Load Sample Info
Sample.Info <- read.csv(file="Respirometry_Data.csv", header=T) #read in volume and sample.info data
Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume

#Merge with sample info
Resp <- merge(Resp.R,Sample.Info, by="Fragment.ID")
Photo <- merge(Photo.R,Sample.Info, by="Fragment.ID")

#Load Surface Area standards from surface area by foil method (Marsh JA (1970) Primary productivity of reef-building calcareous red algae. Ecology 51:255–263)
FoilStnd<-read.csv('foil.standards.csv', header=T, sep=",") #load data
FoilStnd.lm <- lm(FoilStnd$surface.area.cm2 ~ FoilStnd$mass.g) #calculate regression of foil as a function of mass
summary(FoilStnd.lm) #view results
summary(FoilStnd.lm)$r.squared #view R2
FoilStnd.lm$coef #view coefficients

SA.Data <- read.csv('Mcap.surface.area.foil.csv', header=T, sep=",", na.strings="NA") #read in coral foil mass data
SA.Data$SA.final.cm2 <- (FoilStnd.lm$coef[2]*SA.Data$final.mass.g)+FoilStnd.lm$coef[1] #calculate surface area for final live area
SA.Data$SA.dead.cm2 <- (FoilStnd.lm$coef[2]*SA.Data$dead.mass.g)+FoilStnd.lm$coef[1] #calculate surface area for dead area
SA.Data$SA.total.initial.cm2 <- (FoilStnd.lm$coef[2]*(SA.Data$dead.mass.g+SA.Data$final.mass.g))+FoilStnd.lm$coef[1] #calculate total live surface area at the start of the experiment

#Merge with surface area data
Resp <- merge(Resp, SA.Data[,c(1,6)], by="Fragment.ID")
Photo <- merge(Photo, SA.Data[,c(1,6)], by="Fragment.ID")

#Account for blank rate
Amb.P.Blank <- Photo.R[which(Photo.R$Fragment.ID == "AmbientBlank1.P"),3]
Amb.R.Blank <- Resp.R[which(Resp.R$Fragment.ID == "AmbientBlank1.R"),3]
High.P.Blank <- Photo.R[which(Photo.R$Fragment.ID == "CO2Blank1.P"),3]
High.R.Blank <- Resp.R[which(Resp.R$Fragment.ID == "CO2Blank1.R"),3]

#Respiration
Amb.indices <- (Resp$Treatment=="Ambient") # indices is a logical vector with TRUE and FALSE
High.indices <- (Resp$Treatment=="High") # indices is a logical vector with TRUE and FALSE
Resp$corr.rate <- NA
Resp$corr.rate[Amb.indices] <- (Resp$µmol.L.sec - Amb.R.Blank)[Amb.indices]
Resp$corr.rate[High.indices] <- (Resp$µmol.L.sec - High.R.Blank)[High.indices]

#Photosynthesis
Amb.indices <- (Photo$Treatment=="Ambient") # indices is a logical vector with TRUE and FALSE
High.indices <- (Photo$Treatment=="High") # indices is a logical vector with TRUE and FALSE
Photo$corr.rate <- NA
Photo$corr.rate[Amb.indices] <- (Photo$µmol.L.sec - Amb.P.Blank)[Amb.indices]
Photo$corr.rate[High.indices] <- (Photo$µmol.L.sec - High.P.Blank)[High.indices]

#Convert to µmol L-1 min-1
Photo$rate.L.min <- Photo$corr.rate*60*60 #converts to µmol L-1 hr-1
Resp$rate.L.min <- Resp$corr.rate*60*60 #converts to µmol L-1 hr-1

#Convert to µmol min-1
Photo$rate.min <- Photo$rate.L.min*Photo$Vol.L 
Resp$rate.min <- Resp$rate.L.min*Resp$Vol.L 

#Convert to µmol cm-2 min-1 (Normalize to Surface Area)
Photo$rate <- Photo$rate.min/Photo$SA.total.initial.cm2
Resp$rate <- Resp$rate.min/Resp$SA.total.initial.cm2
  


#calculate gross photosynthesis Pnet -- Rdark
resp.data <- merge(Photo[,c(1,4,5,12)],Resp[,c(1,12)],  by="Fragment.ID")
colnames(resp.data)[4] <- "Pnet_umol.cm2.hr"
colnames(resp.data) [5] <- "Rdark_umol.cm2.hr"
resp.data$Pgross_umol.cm2.hr <- resp.data$Pnet_umol.cm2.hr-resp.data$Rdark_umol.cm2.hr

#Calculate and plot by conditions
Pnet.mean <- aggregate(Pnet_umol.cm2.hr ~ History*Treatment, data=resp.data, mean)
Pnet.se <- aggregate(Pnet_umol.cm2.hr ~ History*Treatment, data=resp.data, std.error)
Pnet.n <- aggregate(Pnet_umol.cm2.hr ~ History*Treatment, data=resp.data, length)

Pnet.means <- cbind(Pnet.mean, Pnet.se$Pnet_umol.cm2.hr, Pnet.n$Pnet_umol.cm2.hr) #combine mean and standard error results
colnames(Pnet.means) <- c("History","Treatment", "mean", "se", "n")  #rename columns to describe contents
Pnet.means$Metric <- "Net Photosynthesis"

Rdark.mean <- aggregate(Rdark_umol.cm2.hr ~ History*Treatment, data=resp.data, mean)
Rdark.se <- aggregate(Rdark_umol.cm2.hr ~ History*Treatment, data=resp.data, std.error)
Rdark.n <- aggregate(Rdark_umol.cm2.hr ~ History*Treatment, data=resp.data, length)

Rdark.means <- cbind(Rdark.mean, Rdark.se$Rdark_umol.cm2.hr, Rdark.n$Rdark_umol.cm2.hr) #combine mean and standard error results
colnames(Rdark.means) <- c("History","Treatment", "mean", "se", "n")  #rename columns to describe contents
Rdark.means$Metric <- "Dark Respiration"

Pgross.mean <- aggregate(Pgross_umol.cm2.hr ~ History*Treatment, data=resp.data, mean)
Pgross.se <- aggregate(Pgross_umol.cm2.hr ~ History*Treatment, data=resp.data, std.error)
Pgross.n <- aggregate(Pgross_umol.cm2.hr ~ History*Treatment, data=resp.data, length)

Pgross.means <- cbind(Pgross.mean, Pgross.se$Pgross_umol.cm2.hr, Pgross.n$Pgross_umol.cm2.hr) #combine mean and standard error results
colnames(Pgross.means) <- c("History","Treatment", "mean", "se", "n")  #rename columns to describe contents
Pgross.means$Metric <- "Gross Photosynthesis"

Fig.Pn <-  ggplot(Pnet.means, aes(x=Treatment, y=mean,  group=History)) + #set up plot information
  geom_errorbar(aes(x=Treatment, ymax=mean+se, ymin=mean-se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=History), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #label x axis
  ylab("Net Photosynthesis µmol cm-2 h-1") + #label y axis
  ylim(-1, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position='none') + #set legend location
  ggtitle("Net Photosythesis") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Pn #view plot


Fig.Pg <-  ggplot(Pgross.means, aes(x=Treatment, y=mean,  group=History)) + #set up plot information
  geom_errorbar(aes(x=Treatment, ymax=mean+se, ymin=mean-se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=History), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #label x axis
  ylab("Gross Photosynthesis µmol cm-2 h-1") + #label y axis
  ylim(-1, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position='none') + #set legend location
  ggtitle("Gross Photosythesis") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Pg #view plot

Fig.Rd <-  ggplot(Rdark.means, aes(x=Treatment, y=mean,  group=History)) + #set up plot information
  geom_errorbar(aes(x=Treatment, ymax=mean+se, ymin=mean-se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=History), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #label x axis
  ylab("Dark Respiration µmol cm-2 h-1") + #label y axis
  ylim(-1, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(.5, .7)) + #set legend location
  ggtitle("Dark Respiration") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Rd #view plot

dev.off()
Pn.lm <- lm(Pnet_umol.cm2.hr ~ History*Treatment, data=resp.data)
Pn.res <- anova(Pn.lm)
Pn.res 
par(mfrow=c(1,3))
hist(Pn.lm$residuals)
boxplot(Pn.lm$residuals)
qqplot(Pn.lm$fitted.values, Pn.lm$residuals)

Pg.lm <- lm(Pgross_umol.cm2.hr ~ History*Treatment, data=resp.data)
Pg.res <- anova(Pg.lm)
Pg.res
par(mfrow=c(1,3))
hist(Pg.lm$residuals)
boxplot(Pg.lm$residuals)
qqplot(Pg.lm$fitted.values, Pg.lm$residuals)

Rd.lm <- lm(Rdark_umol.cm2.hr ~ History*Treatment, data=resp.data)
Rd.res <- anova(Rd.lm)
Rd.res
par(mfrow=c(1,3))
hist(Rd.lm$residuals)
boxplot(Rd.lm$residuals)
qqplot(Rd.lm$fitted.values, Rd.lm$residuals)

setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/") #set working
Resp.Figs <- arrangeGrob(Fig.Pn, Fig.Rd, Fig.Pg, ncol=3)
ggsave(file="Respirometry.pdf", Resp.Figs, width = 11, height = 6, units = c("in"))




  
  





