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
file.names<-list.files(path = path.p, pattern = "csv$")
photo.R <- data.frame(matrix(NA, nrow=length(file.names), ncol=3, dimnames=list(file.names,c("Sample.ID","Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data <-read.table(file.path(path.p,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data  <- Photo.Data [,c(9,16)]
  n<-dim(Photo.Data )[1]
  Photo.Data <-Photo.Data [120:(n-3),]
  n<-dim(Photo.Data )[1]
  Photo.Data $sec <- 1:n
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/photo.thinning.pdf")
  par(omi=rep(0.3, 4))
  par(mfrow=c(1,2))
  plot(Value ~ sec, data=Photo.Data , 
       xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $Value ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  
  Photo.Data   <-  thinData(Photo.Data , by=10)$newData1
  Photo.Data $sec <- as.numeric(rownames(Photo.Data ))
  plot(Value ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE)
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
  
  photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts them in the dataframe
  photo.R[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}

photo.R$rate.min <- photo.R$Slope*60
photo.R

#RESPIRATION

path.r<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/Respirometry/Resp/" #the location of all your respirometry files
file.names.r<-list.files(path = path.r, pattern = "csv$")
Resp.R <- data.frame(matrix(NA, nrow=length(file.names.r), ncol=3, dimnames=list(file.names.r,c("Sample.ID","Intercept", "Slope")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names.r)) { # for every file in list start at the first and run this following function
  Resp.Data <-read.table(file.path(path.r,file.names.r[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Resp.Data  <- Resp.Data[,c(9,16)]
  n<-dim(Resp.Data)[1]
  Resp.Data <-Resp.Data[120:(n-3),]
  n<-dim(Resp.Data)[1]
  Resp.Data$sec <- 1:n
  pdf("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/resp.thinning.pdf")
  par(omi=rep(0.3, 4))
  par(mfrow=c(1,2))
  plot(Value ~ sec, data=Resp.Data , 
       xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE)
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
  Resp.R[i,1] <- substr(file.names.r[i],1,8) #stores the file name in the Date column
}

Resp.R$rate.min <- Resp.R$Slope*60
Resp.R


Flux.Info <- read.csv(file="Respirometry_Data.csv", header=T)
Flux.Info$Vol.L <- (Flux.Info$Chamber.Vol.mL-Flux.Info$Sample.Displacement.ml)/1000

# Amb.Photo.Blank
# High.Photo.Blank
# Amb.Resp.Blank 
# High.Resp.Blank



#ANOTHER OPTION
#want to find photosyn local slope
#want to find resp local slope
#want to identify upper and lower data points and 
#then use them to extract all the unthinned data and re-run the regression


#Load calculated data
data <- read.csv("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/Respirometry_Data.csv", header=T, sep=",")

#####RESPIRATION#####
resp.data <- read.csv("Respirometry_Data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

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

resp.data <- merge(resp.data, SA.Data, by="Fragment.ID")

resp.data$Pnet_umol.min <- (resp.data$Pnet_umol.L.min*(resp.data$Chamber.Vol.mL/1000))
resp.data$Rdark_umol.min <- (resp.data$Rdark_umol.L.min*(resp.data$Chamber.Vol.mL/1000))
resp.data$Pnet.Blank_umol.min <- (resp.data$Blank_Pnet_umol.L.m*(resp.data$Chamber.Vol.mL/1000))
resp.data$Rdark.Blank_umol.min <- (resp.data$Blank_Rdark_umol.L.m*(resp.data$Chamber.Vol.mL/1000))
resp.data$Coral_Pnet_umol.min <- resp.data$Pnet_umol.min-resp.data$Pnet.Blank_umol.min
resp.data$Coral_Rdark_umol.min <- resp.data$Rdark_umol.min-resp.data$Rdark.Blank_umol.min
resp.data$Coral_Pnet_umol.min.cm2 <- resp.data$Coral_Pnet_umol.min/resp.data$SA.total.initial.cm2
resp.data$Coral_Rdark_umol.min.cm2 <- resp.data$Coral_Rdark_umol.min/resp.data$SA.total.initial.cm2
resp.data$Coral_Pnet_umol.hr.cm2 <- resp.data$Coral_Pnet_umol.min.cm2*60
resp.data$Coral_Rdark_umol.hr.cm2 <- resp.data$Coral_Rdark_umol.min.cm2*60
resp.data$Coral_Pgross_umol.hr.cm2 <- resp.data$Coral_Pnet_umol.hr.cm2 - resp.data$Coral_Rdark_umol.hr.cm2




Pnet.mean <- aggregate(Coral_Pnet_umol.hr.cm2 ~ History*Treatment, data=resp.data, mean)
Pnet.se <- aggregate(Coral_Pnet_umol.hr.cm2 ~ History*Treatment, data=resp.data, std.error)
Pnet.n <- aggregate(Coral_Pnet_umol.hr.cm2 ~ History*Treatment, data=resp.data, length)

Pnet.means <- cbind(Pnet.mean, Pnet.se$Coral_Pnet_umol.hr.cm2, Pnet.n$Coral_Pnet_umol.hr.cm2) #combine mean and standard error results
colnames(Pnet.means) <- c("History","Treatment", "mean", "se", "n")  #rename columns to describe contents
Pnet.means$Metric <- "Net Photosynthesis"

Rdark.mean <- aggregate(Coral_Rdark_umol.hr.cm2 ~ History*Treatment, data=resp.data, mean)
Rdark.se <- aggregate(Coral_Rdark_umol.hr.cm2 ~ History*Treatment, data=resp.data, std.error)
Rdark.n <- aggregate(Coral_Rdark_umol.hr.cm2 ~ History*Treatment, data=resp.data, length)

Rdark.means <- cbind(Rdark.mean, Rdark.se$Coral_Rdark_umol.hr.cm2, Rdark.n$Coral_Rdark_umol.hr.cm2) #combine mean and standard error results
colnames(Rdark.means) <- c("History","Treatment", "mean", "se", "n")  #rename columns to describe contents
Rdark.means$Metric <- "Dark Respiration"

Pgross.mean <- aggregate(Coral_Pgross_umol.hr.cm2 ~ History*Treatment, data=resp.data, mean)
Pgross.se <- aggregate(Coral_Pgross_umol.hr.cm2 ~ History*Treatment, data=resp.data, std.error)
Pgross.n <- aggregate(Coral_Pgross_umol.hr.cm2 ~ History*Treatment, data=resp.data, length)

Pgross.means <- cbind(Pgross.mean, Pgross.se$Coral_Pgross_umol.hr.cm2, Pgross.n$Coral_Pgross_umol.hr.cm2) #combine mean and standard error results
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
Pn.lm <- lm(Coral_Pnet_umol.hr.cm2 ~ History*Treatment, data=resp.data)
Pn.res <- anova(Pn.lm)
Pn.res 
par(mfrow=c(1,3))
hist(Pn.lm$residuals)
boxplot(Pn.lm$residuals)
qqplot(Pn.lm$fitted.values, Pn.lm$residuals)

Pg.lm <- lm(Coral_Pgross_umol.hr.cm2 ~ History*Treatment, data=resp.data)
Pg.res <- anova(Pg.lm)
Pg.res
par(mfrow=c(1,3))
hist(Pg.lm$residuals)
boxplot(Pg.lm$residuals)
qqplot(Pg.lm$fitted.values, Pg.lm$residuals)

Rd.lm <- lm(Coral_Rdark_umol.hr.cm2 ~ History*Treatment, data=resp.data)
Rd.res <- anova(Rd.lm)
Rd.res
par(mfrow=c(1,3))
hist(Rd.lm$residuals)
boxplot(Rd.lm$residuals)
qqplot(Rd.lm$fitted.values, Rd.lm$residuals)

setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/") #set working
Resp <- arrangeGrob(Fig.Pn, Fig.Rd, Fig.Pg, ncol=3)
ggsave(file="Respirometry.pdf", Resp, width = 11, height = 6, units = c("in"))




  
  





