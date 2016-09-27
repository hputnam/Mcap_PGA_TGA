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
library("MASS")

#Required Data files

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data/") #set working

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

#Fecundity Data
fecundity.data <- read.csv("Mcap_Spawning.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
bun.vol <- 0.0002 #0.0002 ml per bundle based on a sperical bundle with an average diameter of 770 µm from Padilla-Gamino et al 2011
#June Data
fecundity.data$bundle.vol_june1 <- (fecundity.data$Time_20160605_Bun.Num*bun.vol)+fecundity.data$Time_20160605_EggVol
fecundity.data$bundle.vol_june2 <- (fecundity.data$Time_20160606_Bun.Num*bun.vol)+fecundity.data$Time_20160606_EggVol
fecundity.data$bundle.vol_june3 <- (fecundity.data$Time_20160607_Bun.Num*bun.vol)+fecundity.data$Time_20160607_EggVol
fecundity.data$june <- rowSums(fecundity.data[,c("bundle.vol_june1", "bundle.vol_june2", "bundle.vol_june3")], na.rm=T)
fecundity.data$june[fecundity.data$june == 0] <- NA

#July Data
fecundity.data$bundle.vol_july1 <- (fecundity.data$Time_20160704_Bun.Num*bun.vol)+fecundity.data$Time_20160704_EggVol
fecundity.data$bundle.vol_july2 <- (fecundity.data$Time_20160705_Bun.Num*bun.vol)+fecundity.data$Time_20160705_EggVol
fecundity.data$bundle.vol_july3 <- (fecundity.data$Time_20160706_Bun.Num*bun.vol)+fecundity.data$Time_20160706_EggVol
fecundity.data$july <- rowSums(fecundity.data[,c("bundle.vol_july1", "bundle.vol_july2", "bundle.vol_july3")], na.rm=T)
fecundity.data$july[fecundity.data$july == 0] <- NA

#August Data
fecundity.data$bundle.vol_aug1 <- (fecundity.data$Time_20160802_Bun.Num*bun.vol)
fecundity.data$bundle.vol_aug2 <- (fecundity.data$Time_20160803_Bun.Num*bun.vol)
fecundity.data$bundle.vol_aug3 <- (fecundity.data$Time_20160804_Bun.Num*bun.vol)
fecundity.data$aug <- rowSums(fecundity.data[,c("bundle.vol_aug1", "bundle.vol_aug2", "bundle.vol_aug3")], na.rm=T)
fecundity.data$aug[fecundity.data$aug == 0] <- NA

fecund <- merge(fecundity.data, SA.Data)
fecund$june.ml.SA <- (fecund$june/fecund$SA.total.initial.cm2)*10^4
fecund$july.ml.SA <- (fecund$july/fecund$SA.total.initial.cm2)*10^4
fecund$aug.ml.SA <- (fecund$aug/fecund$SA.total.initial.cm2)*10^4
fecund$tot.ml.SA <- rowSums(fecund[,c("june.ml.SA", "july.ml.SA", "aug.ml.SA")], na.rm=T)
fecund$tot.ml.SA[fecund$tot.ml.SA == 0.000] <- NA
fecund$spawn[fecund$tot.ml.SA > 0.000] <- 1
fecund$spawn[is.na(fecund$tot.ml.SA)] <- 0


#June
june.fec.mean <- aggregate(june.ml.SA ~ Parental.Performance*Treatment, data=fecund, mean)
june.fec.se <- aggregate(june.ml.SA ~ Parental.Performance*Treatment, data=fecund, std.error)
june.fec.n <- aggregate(june.ml.SA ~ Parental.Performance*Treatment, data=fecund, length)

june.fec.means <- cbind(june.fec.mean, june.fec.se$june.ml.SA, june.fec.n$june.ml.SA)
colnames(june.fec.means) <-c("Parental.Performance",  "Treatment",	"mean", "se", "n")
june.fec.means$TP <- paste(june.fec.means$Parental.Performance, june.fec.means$Treatment, sep='') #add a concatenated history by treatment grouping

Fig.June.Fecundity <- ggplot(june.fec.means, aes(x=Treatment, y=mean, group=Parental.Performance)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #Label the X Axis
  ylab("June Fecundity ml cm-2 10^4") + #Label the Y Axis
  #ylim(0.6, 0.75) + #set Y limits
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig.June.Fecundity 

#july
july.fec.mean <- aggregate(july.ml.SA ~ Parental.Performance*Treatment, data=fecund, mean)
july.fec.se <- aggregate(july.ml.SA ~ Parental.Performance*Treatment, data=fecund, std.error)
july.fec.n <- aggregate(july.ml.SA ~ Parental.Performance*Treatment, data=fecund, length)

july.fec.means <- cbind(july.fec.mean, july.fec.se$july.ml.SA, july.fec.n$july.ml.SA)
colnames(july.fec.means) <-c("Parental.Performance",  "Treatment",  "mean", "se", "n")
july.fec.means$TP <- paste(july.fec.means$Parental.Performance, july.fec.means$Treatment, sep='') #add a concatenated history by treatment grouping

Fig.july.Fecundity <- ggplot(july.fec.means, aes(x=Treatment, y=mean, group=Parental.Performance)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #Label the X Axis
  ylab("July Fecundity ml cm-2 10^4") + #Label the Y Axis
  #ylim(0.6, 0.75) + #set Y limits
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig.july.Fecundity 

#aug
aug.fec.mean <- aggregate(aug.ml.SA ~ Parental.Performance*Treatment, data=fecund, mean)
aug.fec.se <- aggregate(aug.ml.SA ~ Parental.Performance*Treatment, data=fecund, std.error)
aug.fec.n <- aggregate(aug.ml.SA ~ Parental.Performance*Treatment, data=fecund, length)

aug.fec.means <- cbind(aug.fec.mean, aug.fec.se$aug.ml.SA, aug.fec.n$aug.ml.SA)
colnames(aug.fec.means) <-c("Parental.Performance",  "Treatment",  "mean", "se", "n")
aug.fec.means$TP <- paste(aug.fec.means$Parental.Performance, aug.fec.means$Treatment, sep='') #add a concatenated history by treatment grouping

Fig.aug.Fecundity <- ggplot(aug.fec.means, aes(x=Treatment, y=mean, group=Parental.Performance)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #Label the X Axis
  ylab("aug Fecundity ml cm-2 10^4") + #Label the Y Axis
  #ylim(0.6, 0.75) + #set Y limits
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig.aug.Fecundity 

#tot
tot.fec.mean <- aggregate(tot.ml.SA ~ Parental.Performance*Treatment, data=fecund, mean)
tot.fec.se <- aggregate(tot.ml.SA ~ Parental.Performance*Treatment, data=fecund, std.error)
tot.fec.n <- aggregate(tot.ml.SA ~ Parental.Performance*Treatment, data=fecund, length)

tot.fec.means <- cbind(tot.fec.mean, tot.fec.se$tot.ml.SA, tot.fec.n$tot.ml.SA)
colnames(tot.fec.means) <-c("Parental.Performance",  "Treatment",  "mean", "se", "n")
tot.fec.means$TP <- paste(tot.fec.means$Parental.Performance, tot.fec.means$Treatment, sep='') #add a concatenated history by treatment grouping

Fig.tot.Fecundity <- ggplot(tot.fec.means, aes(x=Treatment, y=mean, group=Parental.Performance)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Treatment") + #Label the X Axis
  ylab(expression(Fecundity~ml~cm^{-2}~x~10^{4})) + #Label the Y Axis
  #expression(Concentration~mg~L^{-1}))
  ylim(0, 2) + #set Y limits
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(.5, .8)) + #set legend location 
  ggtitle("B) Fecundity") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                size = 12, 
                                hjust = 0)) #set title attributes
Fig.tot.Fecundity 

dev.off()
Fec.lm <- lm(log10(tot.ml.SA) ~ Parental.Performance*Treatment, data=fecund)
Fec.res <- anova(Fec.lm)
Fec.res 
par(mfrow=c(1,3))
hist(Fec.lm$residuals)
boxplot(Fec.lm$residuals)
qqplot(Fec.lm$fitted.values, Fec.lm$residuals)

fecund$tot.ml.SA

tbl <- table(fecund$spawn, fecund$Treatment, fecund$Parental.Performance) #generate a spawning table
tbl #view table
chisq.test(tbl) #chi square test of independence H0 there is no difference in spawning (1) and nonspawning (0) between bleached and unbleached histories and ambient and high pCO2


spawning.freq <- aggregate(spawn ~ Treatment*Parental.Performance, data=fecund, sum)
spawning.count <- aggregate(spawn ~ Treatment*Parental.Performance, data=fecund, length)
spawning <-cbind(spawning.freq, spawning.count$spawn)
colnames(spawning) <- c("Treatment", "Parental.Performance", "count", "total")
spawning$prop <-spawning$count/spawning$total
spawning$Spawn <- "Yes"

no.spawn <-(1-spawning$prop)*spawning$total
no.spawn <- data.frame(spawning.freq$Treatment, spawning.freq$Parental.Performance,no.spawn, spawning$total)
no.spawn
colnames(no.spawn) <- c("Treatment", "Parental.Performance", "count", "total")
no.spawn$prop <-no.spawn$count/no.spawn$total
no.spawn$Spawn <- "No"

spawning.data <- rbind(spawning, no.spawn)
spawning.data$TP <- paste(spawning.data$Parental.Performance, spawning.data$Treatment, sep='_') #add a concatenated history by treatment grouping

Fig.Spawning <- ggplot(spawning.data, aes(x=TP, y=prop, fill=Spawn)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("gray", "black")) +
  xlab("Treatment") + #Label the X Axis
  ylab("Proportion Spawning") + #Label the Y Axis
  #scale_x_discrete(labels=c(Ambient pH,)) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank(),  #Set plot legend key
        legend.position="top") + #set legend location
  ggtitle("A) Proportion Spawning") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Spawning


##### Larval Survivorship #####
##### Larval Counts #####
larval.counts <- read.csv("Larval_Survivorship.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
larval.counts$Avg.Count <- rowMeans(larval.counts[,c("Count1", "Count2", "Count3")], na.rm = TRUE) #calculate average of counts
larval.counts$larvae.per.ml <- larval.counts$Avg.Count/larval.counts$Sample.vol #calculate density
larval.counts$total.larvae <- larval.counts$Avg.Count/(larval.counts$Sample.vol/larval.counts$Tripour.vol)

#Survivorship
Survive <- reshape(larval.counts, idvar = "Conical", timevar = "Time.Point", direction = "wide")
Survive$Percent.Sur <- 100- ((Survive$total.larvae.Time1 - Survive$total.larvae.Time2)/Survive$total.larvae.Time1)*100

#Treatments
gmean_Sur <- tapply(Survive$Percent.Sur, list(Survive$Group.Time1), mean, na.rm = TRUE)
gse_Sur <- tapply(Survive$Percent.Sur , list(Survive$Group.Time1), std.error, na.rm = TRUE)
gmean_Sur <- as.data.frame(cbind(gmean_Sur, gse_Sur))

sur.mean <- aggregate(Percent.Sur ~ Parental.Performance.Time1*Parent.Treatment.Time1*Larval.Treatment.Time1, data=Survive, mean)
sur.se <- aggregate(Percent.Sur ~ Parental.Performance.Time1*Parent.Treatment.Time1*Larval.Treatment.Time1, data=Survive, std.error)
sur.n <- aggregate(Percent.Sur ~ Parental.Performance.Time1*Parent.Treatment.Time1*Larval.Treatment.Time1, data=Survive, length)

Sur.means <- cbind(sur.mean, sur.se$Percent.Sur, sur.n$Percent.Sur)
colnames(Sur.means) <-c("Parental.Performance",  "Parent.Treatment",  "Larval.Treatment",	"mean", "se", "n")


Fig.Sur <-  ggplot(Sur.means, aes(x=Larval.Treatment, y=mean,  group=Parental.Performance)) + #set up plot information
  geom_errorbar(aes(x=Larval.Treatment, ymax=mean+se, ymin=mean-se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Larval Treatment") + #label x axis
  ylab("% Survivorship") + #label y axis
  ylim(0,35)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position='none') + #set legend location 
  ggtitle("C) Larval Survivorship") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
Fig.Sur #view plot

Sur.lm <- lm(sqrt(Percent.Sur) ~ Parental.Performance.Time1*Larval.Treatment.Time1, data=Survive)
Sur.res <- anova(Sur.lm)
Sur.res
par(mfrow=c(1,3))
hist(Sur.lm$residuals)
boxplot(Sur.lm$residuals)
qqplot(Sur.lm$fitted.values, Sur.lm$residuals)



setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/") #set working
Repro <- arrangeGrob(Fig.Spawning, Fig.tot.Fecundity, Fig.Sur, ncol=3)
ggsave(file="Reproduction.pdf", Repro, width = 11, height = 6, units = c("in"))


