#Title: PGA TGA PAM Dark Adapted Yield Data
#Date Last Modified: 20160331
#See Readme file for details

rm(list=ls()) #clears workspace 


#Read in required libraries
##### Include Versions of libraries
library('plyr')
library('ggplot2')
library("sciplot")
library("plotrix")
library("reshape2")
library("nlme") #mixed model, repeated measures ANOVA
library("lsmeans") #mixed model posthoc  statistical comparisons

#Required Data files
#PGA_TGA_PAM2016.csv

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data") #set working

#Load data
Time0 <- read.csv("Time0_20160228_PAM.csv", header=T) #load data
Time1 <- read.csv("Time1_20160310_PAM.csv", header=T) #load data
Time2 <- read.csv("Time2_20160317_PAM.csv", header=T) #load data
Time3 <- read.csv("Time3_20160324_PAM.csv", header=T) #load data
Time4 <- read.csv("Time4_20160407_PAM.csv", header=T) #load data
Time5 <- read.csv("Time5_20160421_PAM.csv", header=T) #load data
Time6 <- read.csv("Time6_20160505_PAM.csv", header=T) #load data
Time7 <- read.csv("Time7_20160519_PAM.csv", header=T) #load data
Time8 <- read.csv("Time8_20160616_PAM.csv", header=T) #load data
Time9 <- read.csv("Time9_20160714_PAM.csv", header=T) #load data
Time10 <- read.csv("Time10_20160811_PAM.csv", header=T) #load data


#check yield calculations versus records and build dataframe
checks <- (Time0[,1:5]) #add sample info to dataframe
checks$T0.Yield.Check.Rep1 <- Time0$Yield.Rep1 - ((Time0$Fm.Rep1-Time0$Fo.Rep1)/Time0$Fm.Rep1)
checks$T0.Yield.Check.Rep2 <- Time0$Yield.Rep2 - ((Time0$Fm.Rep2-Time0$Fo.Rep2)/Time0$Fm.Rep2)
checks$T1.Yield.Check.Rep1 <- Time1$Yield.Rep1 - ((Time1$Fm.Rep1-Time1$Fo.Rep1)/Time1$Fm.Rep1)
checks$T1.Yield.Check.Rep2 <- Time1$Yield.Rep2 - ((Time1$Fm.Rep2-Time1$Fo.Rep2)/Time1$Fm.Rep2)
checks$T2.Yield.Check.Rep1 <- Time2$Yield.Rep1 - ((Time2$Fm.Rep1-Time2$Fo.Rep1)/Time2$Fm.Rep1)
checks$T2.Yield.Check.Rep2 <- Time2$Yield.Rep2 - ((Time2$Fm.Rep2-Time2$Fo.Rep2)/Time2$Fm.Rep2)
checks$T3.Yield.Check.Rep1 <- Time3$Yield.Rep1 - ((Time3$Fm.Rep1-Time3$Fo.Rep1)/Time3$Fm.Rep1)
checks$T3.Yield.Check.Rep2 <- Time3$Yield.Rep2 - ((Time3$Fm.Rep2-Time3$Fo.Rep2)/Time3$Fm.Rep2)
checks$T4.Yield.Check.Rep1 <- Time4$Yield.Rep1 - ((Time4$Fm.Rep1-Time4$Fo.Rep1)/Time4$Fm.Rep1)
checks$T4.Yield.Check.Rep2 <- Time4$Yield.Rep2 - ((Time4$Fm.Rep2-Time4$Fo.Rep2)/Time4$Fm.Rep2)
checks$T5.Yield.Check.Rep1 <- Time5$Yield.Rep1 - ((Time5$Fm.Rep1-Time5$Fo.Rep1)/Time5$Fm.Rep1)
checks$T5.Yield.Check.Rep2 <- Time5$Yield.Rep2 - ((Time5$Fm.Rep2-Time5$Fo.Rep2)/Time5$Fm.Rep2)
checks$T6.Yield.Check.Rep1 <- Time6$Yield.Rep1 - ((Time6$Fm.Rep1-Time6$Fo.Rep1)/Time6$Fm.Rep1)
checks$T6.Yield.Check.Rep2 <- Time6$Yield.Rep2 - ((Time6$Fm.Rep2-Time6$Fo.Rep2)/Time6$Fm.Rep2)
checks$T7.Yield.Check.Rep1 <- Time7$Yield.Rep1 - ((Time7$Fm.Rep1-Time7$Fo.Rep1)/Time7$Fm.Rep1)
checks$T7.Yield.Check.Rep2 <- Time7$Yield.Rep2 - ((Time7$Fm.Rep2-Time7$Fo.Rep2)/Time7$Fm.Rep2)
checks$T8.Yield.Check.Rep1 <- Time8$Yield.Rep1 - ((Time8$Fm.Rep1-Time8$Fo.Rep1)/Time8$Fm.Rep1)
checks$T8.Yield.Check.Rep2 <- Time8$Yield.Rep2 - ((Time8$Fm.Rep2-Time8$Fo.Rep2)/Time8$Fm.Rep2)
checks$T9.Yield.Check.Rep1 <- Time9$Yield.Rep1 - ((Time9$Fm.Rep1-Time9$Fo.Rep1)/Time9$Fm.Rep1)
checks$T9.Yield.Check.Rep2 <- Time9$Yield.Rep2 - ((Time9$Fm.Rep2-Time9$Fo.Rep2)/Time9$Fm.Rep2)
checks$T10.Yield.Check.Rep1 <- Time10$Yield.Rep1 - ((Time10$Fm.Rep1-Time10$Fo.Rep1)/Time10$Fm.Rep1)
checks$T10.Yield.Check.Rep2 <- Time10$Yield.Rep2 - ((Time10$Fm.Rep2-Time10$Fo.Rep2)/Time10$Fm.Rep2)

#Rep1
par(mfrow=c(3,4), oma=c(3,0,0,0))
boxplot(checks$T0.Yield.Check.Rep1, main="Time0") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T1.Yield.Check.Rep1, main="Time1") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T2.Yield.Check.Rep1, main="Time2") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T3.Yield.Check.Rep1, main="Time3") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T4.Yield.Check.Rep1, main="Time4") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T5.Yield.Check.Rep1, main="Time5") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T6.Yield.Check.Rep1, main="Time6") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T7.Yield.Check.Rep1, main="Time7") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T8.Yield.Check.Rep1, main="Time8") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T9.Yield.Check.Rep1, main="Time9") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T10.Yield.Check.Rep1, main="Time10") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
mtext("Yield Calculation Check Rep1", outer=T, cex = 1.5, side=1)

#list which calculations are >0.002 from written value
which(abs(checks$T0.Yield.Check.Rep1) > 0.002)
which(abs(checks$T1.Yield.Check.Rep1) > 0.002)
which(abs(checks$T2.Yield.Check.Rep1) > 0.002)
which(abs(checks$T3.Yield.Check.Rep1) > 0.002)
which(abs(checks$T4.Yield.Check.Rep1) > 0.002)
which(abs(checks$T5.Yield.Check.Rep1) > 0.002)
which(abs(checks$T6.Yield.Check.Rep1) > 0.002)
which(abs(checks$T7.Yield.Check.Rep1) > 0.002)
which(abs(checks$T8.Yield.Check.Rep1) > 0.002)
which(abs(checks$T9.Yield.Check.Rep1) > 0.002)
which(abs(checks$T10.Yield.Check.Rep1) > 0.002)

#Rep2
par(mfrow=c(3,4), oma=c(3,0,0,0))
boxplot(checks$T0.Yield.Check.Rep2, main="Time0") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T1.Yield.Check.Rep2, main="Time1") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T2.Yield.Check.Rep2, main="Time2") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T3.Yield.Check.Rep2, main="Time3") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T4.Yield.Check.Rep2, main="Time4") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T5.Yield.Check.Rep2, main="Time5") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T6.Yield.Check.Rep2, main="Time6") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T7.Yield.Check.Rep2, main="Time7") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T8.Yield.Check.Rep2, main="Time8") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T9.Yield.Check.Rep2, main="Time9") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
boxplot(checks$T10.Yield.Check.Rep2, main="Time10") #check for outliers
abline(h=c(0.002,-0.002), col = "red")
mtext("Yield Calculation Check Rep2", outer=T, cex = 1.5, side=1)

#list which calculations are >0.002 from written value
which(abs(checks$T0.Yield.Check.Rep2) > 0.002)
which(abs(checks$T1.Yield.Check.Rep2) > 0.002)
which(abs(checks$T2.Yield.Check.Rep2) > 0.002)
which(abs(checks$T3.Yield.Check.Rep2) > 0.002)
which(abs(checks$T4.Yield.Check.Rep2) > 0.002)
which(abs(checks$T5.Yield.Check.Rep2) > 0.002)
which(abs(checks$T6.Yield.Check.Rep2) > 0.002)
which(abs(checks$T7.Yield.Check.Rep2) > 0.002)
which(abs(checks$T8.Yield.Check.Rep2) > 0.002)
which(abs(checks$T9.Yield.Check.Rep2) > 0.002)
which(abs(checks$T10.Yield.Check.Rep2) > 0.002)

#Examine rep1 in relation to rep2
par(mfrow=c(3,4), oma=c(3,0,0,0))
plot(Time0$Yield.Rep1,Time0$Yield.Rep2, main="Time0")
T0 <- cor(Time0[,c(10,13)], method = c("pearson"), use = "complete")
T0 <- "NA"
plot(Time1$Yield.Rep1,Time1$Yield.Rep2, main="Time1")
T1 <- cor(Time1[,c(10,13)], method = c("pearson"), use = "complete")
T1 <- T1[2,1]
plot(Time2$Yield.Rep1,Time2$Yield.Rep2, main="Time2")
T2 <- cor(Time2[,c(10,13)], method = c("pearson"), use = "complete")
T2 <- T2[2,1]
plot(Time3$Yield.Rep1,Time3$Yield.Rep2, main="Time3")
T3 <- cor(Time3[,c(10,13)], method = c("pearson"), use = "complete")
T3 <- T3[2,1]
plot(Time4$Yield.Rep1,Time4$Yield.Rep2, main="Time4")
T4 <- cor(Time4[,c(10,13)], method = c("pearson"), use = "complete")
T4 <- T4[2,1]
plot(Time5$Yield.Rep1,Time5$Yield.Rep2, main="Time5")
T5 <- cor(Time5[,c(10,13)], method = c("pearson"), use = "complete")
T5 <- "NA"
plot(Time6$Yield.Rep1,Time6$Yield.Rep2, main="Time6")
T6 <- cor(Time6[,c(10,13)], method = c("pearson"), use = "complete")
T6 <- T6[2,1]
plot(Time7$Yield.Rep1,Time7$Yield.Rep2, main="Time7")
T7 <- cor(Time7[,c(10,13)], method = c("pearson"), use = "complete")
T7 <- T7[2,1]
plot(Time8$Yield.Rep1,Time8$Yield.Rep2, main="Time8")
T8 <- cor(Time8[,c(10,13)], method = c("pearson"), use = "complete")
T8 <- T8[2,1]
plot(Time9$Yield.Rep1,Time9$Yield.Rep2, main="Time9")
T9 <- cor(Time9[,c(10,13)], method = c("pearson"), use = "complete")
T9 <- T9[2,1]
plot(Time10$Yield.Rep1,Time10$Yield.Rep2, main="Time10")
T10 <- cor(Time10[,c(10,13)], method = c("pearson"), use = "complete")
T10 <- T10[2,1]
mtext("Correlation between Rep1 and Rep2", outer=T, cex = 1.5, side=1)
corrs <- rbind(T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10)

#replace outliers with NA
pam.mn <- 0.500
Time0$Yield.Rep1[Time0$Yield.Rep1 < pam.mn] <- NA
Time1$Yield.Rep1[Time1$Yield.Rep1 < pam.mn] <- NA
Time2$Yield.Rep1[Time2$Yield.Rep1 < pam.mn] <- NA
Time3$Yield.Rep1[Time3$Yield.Rep1 < pam.mn] <- NA
Time4$Yield.Rep1[Time4$Yield.Rep1 < pam.mn] <- NA
Time5$Yield.Rep1[Time5$Yield.Rep1 < pam.mn] <- NA
Time6$Yield.Rep1[Time6$Yield.Rep1 < pam.mn] <- NA
Time7$Yield.Rep1[Time7$Yield.Rep1 < pam.mn] <- NA
Time8$Yield.Rep1[Time8$Yield.Rep1 < pam.mn] <- NA
Time9$Yield.Rep1[Time9$Yield.Rep1 < pam.mn] <- NA
Time10$Yield.Rep1[Time10$Yield.Rep1 < pam.mn] <- NA

Time0$Yield.Rep2[Time0$Yield.Rep2 < pam.mn] <- NA
Time1$Yield.Rep2[Time1$Yield.Rep2 < pam.mn] <- NA
Time2$Yield.Rep2[Time2$Yield.Rep2 < pam.mn] <- NA
Time3$Yield.Rep2[Time3$Yield.Rep2 < pam.mn] <- NA
Time4$Yield.Rep2[Time4$Yield.Rep2 < pam.mn] <- NA
Time5$Yield.Rep2[Time5$Yield.Rep2 < pam.mn] <- NA
Time6$Yield.Rep2[Time6$Yield.Rep2 < pam.mn] <- NA
Time7$Yield.Rep2[Time7$Yield.Rep2 < pam.mn] <- NA
Time8$Yield.Rep2[Time8$Yield.Rep2 < pam.mn] <- NA
Time9$Yield.Rep2[Time9$Yield.Rep2 < pam.mn] <- NA
Time10$Yield.Rep2[Time10$Yield.Rep2 < pam.mn] <- NA

#calculate averages and build dataframe
PAM.data <- (Time0[,1:5])
PAM.data$T0.yield <- rowMeans(Time0[,c(10,13)], na.rm = TRUE)
PAM.data$T1.yield <- rowMeans(Time1[,c(10,13)], na.rm = TRUE)
PAM.data$T2.yield <- rowMeans(Time2[,c(10,13)], na.rm = TRUE)
PAM.data$T3.yield <- rowMeans(Time3[,c(10,13)], na.rm = TRUE)
PAM.data$T4.yield <- rowMeans(Time4[,c(10,13)], na.rm = TRUE)
PAM.data$T5.yield <- rowMeans(Time5[,c(10,13)], na.rm = TRUE)
PAM.data$T6.yield <- rowMeans(Time6[,c(10,13)], na.rm = TRUE)
PAM.data$T7.yield <- rowMeans(Time7[,c(10,13)], na.rm = TRUE)
PAM.data$T8.yield <- rowMeans(Time8[,c(10,13)], na.rm = TRUE)
PAM.data$T9.yield <- rowMeans(Time9[,c(10,13)], na.rm = TRUE)
PAM.data$T10.yield <- rowMeans(Time10[,c(10,13)], na.rm = TRUE)
PAM.data[ is.na(PAM.data) ] <- NA

#Check for outliers in averages
mins <- sapply(na.omit(PAM.data[,6:16]), FUN=min)
maxes <- sapply(na.omit(PAM.data[,6:16]), FUN=max)
ranges <- rbind(mins, maxes)
ranges


#view data for outliers
dev.off()
par(mfrow=c(4,3))
boxplot(PAM.data$T0.yield)
boxplot(PAM.data$T1.yield)
boxplot(PAM.data$T2.yield)
boxplot(PAM.data$T3.yield)
boxplot(PAM.data$T4.yield)
boxplot(PAM.data$T5.yield)
boxplot(PAM.data$T6.yield)
boxplot(PAM.data$T7.yield)
boxplot(PAM.data$T8.yield)
boxplot(PAM.data$T9.yield)
boxplot(PAM.data$T10.yield)

#Reshape data
PAM.long <- melt(PAM.data) #reshape data
colnames(PAM.long) <- c("Parental.Performance",  "Parent.ID",  "Fragment.ID",	"Cohort",	"Treatment",	"Time",	"Yield") 
PAM.long <- na.omit(PAM.long) #remove na

#Calculate Summary Statistics
Yield <-ddply(PAM.long, .(Parental.Performance, Treatment, Time), summarize, #For each subset of a data frame, apply function then combine results into a data frame.
              N=length(na.omit(Yield)), # the sample size of the light column summarized by Tank number
              mean = (mean(Yield)),       #take the average of the light column summarized by Tank number
              sd = sd(Yield), # the standard deviation of the temp column summarized by Tank number
              sem = (sd(Yield)/sqrt(N))) #calculate the SEM as the sd/sqrt of the count or data length
Yield$TP <- paste(Yield$Parental.Performance, Yield$Treatment, sep='') #add a concatenated history by treatment grouping
Yield$Date <- c("Feb28", "March10", "March17", "March24", "April07", "April21", "May05", "May19", "June16", "July14", "August11", "Feb28", "March10", "March17", "March24", "April07", "April21", "May05", "May19", "June16", "July14", "August11", "Feb28", "March10", "March17", "March24", "April07", "April21", "May05", "May19", "June16", "July14", "August11", "Feb28", "March10", "March17", "March24", "April07", "April21", "May05", "May19", "June16", "July14", "August11")



#plot the yeild over time with pH indicated lines and parental history by symbols
Fig.PAM <- ggplot(Yield, aes(x=Time, y=mean, group=TP)) + 
  geom_errorbar(aes(ymin=Yield$mean-Yield$sem, ymax=Yield$mean+Yield$sem), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Fv/Fm") + #Label the Y Axis
  ylim(0.6, 0.75) + #set Y limits
  scale_x_discrete(labels=c("Feb28", "March10", "March17", "March24", "April07", "April21", "May05", "May19", "June16", "July14", "August11")) +
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig.PAM


setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/") #set working
ggsave(file="PAM.pdf", Fig.PAM, width = 11, height = 6, units = c("in"))




###Repeated Measures ANOVA
PAM.RM.lme <- lme(Yield ~ Time*Parental.Performance*Treatment, random = ~ Time|Fragment.ID, data=na.omit(PAM.long)) #repeated measures ANOVA with random intercept but not slope (clonal fragments expect to respond the same)
PAM.results <- summary(PAM.RM.lme) #view RM ANOVA summary
PAM.stats <- anova(PAM.RM.lme) #view F and p values
PAM.stats

PAM.RM.posthoc <- lsmeans(PAM.RM.lme, specs=c("variable","Parental.Performance")) #calculate MS means
PAM.RM.posthoc #view results
PAM.RM.posthoc.p <- contrast(PAM.RM.posthoc, method="pairwise", by=c("variable")) #contrast treatment groups within a species at each time point
PAM.RM.posthoc.p #view results

###Testing ANOVA Assumptions
###Data are normally distributed
hist(PAM.RM.lme$residuals) #histogram
qqnorm(PAM.RM.lme$residuals) # normal quantile plot


