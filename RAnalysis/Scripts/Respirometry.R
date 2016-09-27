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

#Required Data files

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/") #set working
mainDir<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/" #set main directory
path<-"/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/Respirometry" #the location of all your respirometry files

# #find all the data files
# file.names<-list.files(path=path) #list all the file names in your data and sample directory
# file.names <- file.names[grep("[.]csv", file.names)] # select only the csv files
# 
# #create an empty dataframe to put the Slope values in
# nrow<-length(file.names) #set number of rows to the number of samples
# Slopes <- matrix(nrow = nrow, ncol = 4) #set the dimensions of the dataframe
# rownames(Slopes)<-file.names #identify row names
# colnames(Slopes)<-c('Sample.ID','Photo','Resp','Eqs') #identify column names
# 
# #run a for loop to bring in the respirometry files one at a time and calculate 2 slopes per file
# for(i in 1: length(file.names)) {
#   Data<-read.csv(text=paste0(head(readLines(file.names[i]), -4), collapse="\n"), header=T, sep=",", na.string="NA",as.is=T)
#   Data<-Data[,c(1,2,7,9,16,20)] #remove unnecessary columns
#   
#   #name of the file without .csv
#   name<-unlist(strsplit(file.names[i], split='.', fixed=TRUE))[1]
#   
#   #find the max to indicate break point
#   #Data$delta_t[which.max(Data$Value)]
#   #identify the linear model
#   lin.mod <- lm(Value~delta_t, data=Data)
# 
#   #Calculate Slopes
#   Slopes[i,1]<-name #add sample name to data output
#   Slopes[i,2]<-lin.mod$coefficients[1] #add y intercept
#   Slopes[i,3]<-lin.mod$coefficients[3] #add slope
# 
# }
# 
# #output in µmol per liter oxygen per min
# Slopes <- data.frame(Slopes)
# 
# class(Slopes)
# 
# par(mfrow=c(6,6))
# par(mar=c(0,1,1,1))
# for(i in 1: length(file.names)) {
#   Data<-read.csv(text=paste0(head(readLines(file.names[i]), -4), collapse="\n"), header=T, sep=",", na.string="NA",as.is=T)
#   plot(Value~delta_t, pch=16, cex=0.5, data=Data)
# }
# 
# 
# par(mfrow=c(6,6))
# par(mar=c(0,1,1,1))
# for(i in 1: length(file.names)) {
#   Data<-read.csv(text=paste0(head(readLines(file.names[i]), -4), collapse="\n"), header=T, sep=",", na.string="NA",as.is=T)
#   plot(Value~delta_t, pch=16, cex=0.5, data=Data)
#   plot(Slopes, add=T, color="red")
# }
# 
# #Load sample data including chamber volume and sample displacement
# #remove data and underscore from ID
# #merge with sample ID
# 
# #Need to subtract blanks
# 
# #Need to account for volume displacement in rate
# 
# #Need to normalize to surface area

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




  
  





