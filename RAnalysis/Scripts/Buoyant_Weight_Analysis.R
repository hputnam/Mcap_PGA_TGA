#Title: PGA TGA Buoyant Weight Data
#Author: HM Putnam
#Date Last Modified: 20160926
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
#Bouyant_weight_calibration_curve.csv
#foil.standards.csv
#Mcap.surface.area.foil.csv

# Set Working Directory:
setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Data") #set working

# load data 
CalData<-read.csv('Bouyant_weight_calibration_curve.csv', header=T, sep=",")
#Data with Standard weights in grams in deionized water, saltwater, and the temp that they were measured at

#plot relationship between temp and Standard in fresh and salt water
plot(CalData$Temp, CalData$StandardSalt, col = 'red', ylab = 'Salt Weight (g)', xlab = 'Temp C')
par(new=TRUE)
plot(CalData$Temp, CalData$StandardFresh, col = 'blue',xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("'Fresh Weight (g)'",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("Salt","Fresh"), bty = "n")

#linear regression between temp and Standard
#Salt water standard measures
StandardSaltModel <- lm(CalData$StandardSalt~CalData$Temp)
summary(StandardSaltModel)
summary(StandardSaltModel)$r.squared
StandardSaltModel$coef

#Fresh water standard measures
StandardFreshModel <- lm(CalData$StandardFresh~CalData$Temp)
summary(StandardFreshModel)
summary(StandardFreshModel)$r.squared
StandardFreshModel$coef

#load coral weight data
BW.data <-read.csv('Mcap_Buoyant_Weight.csv', header=T, na.strings = "NA") 
BW <- BW.data$Mass.g #assign mass data to BW for function
Temp <- BW.data$Temp.C #assign temp to temp for function

BWCalc <- function(StandardAir= 39.108, Temp, BW, CoralDensity = 2.03){
  #Parameters---------------------------------------------------------------
  # StandardAir is the weight of the Standard in air (Default set at 39.108 grams weighed on balance)
  # Temp is the temperature of the water during the coral measurement
  # BW is the buoywant weight of the coral
  # CoralDensity <- 2.03 #g cm-3, set aragonite density for Montipora from literature 
  #Montipora Arag = 2.03 g cm-3 average from table 2 in Anthony 2003 Functional Ecology 17:246-259
  # You can change the density to literatrue values for the specific coral species of interest in the function. 
  
  #Calculation------------------------------------------------------------
  #Step 1: correct the Standard weights for temperature
  # StandardFresh is the weight of the Standard in  fresh water at temp x
  # StandardSalt is the weight of the Standard in salt water at temp x
  
  # This is based on a calibration curve for Standards weighed in fresh and salt water at many temps
StandardFresh <- StandardFreshModel$coef[-1]*Temp + StandardFreshModel$coef[1] 
StandardSalt <- StandardSaltModel$coef[-1]*Temp + StandardSaltModel$coef[1] 
  
  # Step 2: Use weight in air and freshwater of the glass Standard to calculate
  # the density of the Standard (Spencer Davies eq. 4)
FreshDensity <- 1 #Fresh water has a density of 1 g/cm3
StandardDensity <- FreshDensity/(1-(StandardFresh/StandardAir)) 
  
  # Step 3: Calculate the density of seawater using the density of the Standard
  # (Spencer Davies eq. 3)
SWDensity <- StandardDensity*(1-(StandardSalt/StandardAir))
  
  # Step 4: Calculate the dry weight of the coral (Spencer Davies eq. 1)
CoralWeight <- BW/(1-(SWDensity/CoralDensity))
  
  return(CoralWeight) #returns coral dry weights in g
}
  
BW.data$Dry.Weigh.g <- BWCalc(BW=BW, Temp=Temp) #use function to calculate dry weight
hist(BW.data$Dry.Weigh.g) #examine data
boxplot(BW.data$Dry.Weigh.g) #examine data

data <- reshape(BW.data, timevar = "TimePoint", drop = c("Tank", "QC", "Mass.g","Temp.C"), idvar=c("Fragment.ID","Parental.Performance", "Treatment", "Parent.ID","Cohort" ), direction="wide")

#####Growth Normalized to timepoint minus 1 mass#####
# data$time1.growth <- (((data$Dry.Weigh.g.Time1 - data$Dry.Weigh.g.Time0)/data$Dry.Weigh.g.Time0)*100)/data$Days.Time1 #calculate growth rate for time1
# data$time2.growth <- (((data$Dry.Weigh.g.Time2 - data$Dry.Weigh.g.Time1)/data$Dry.Weigh.g.Time1)*100)/data$Days.Time2 #calculate growth rate for time2
# data$time3.growth <- (((data$Dry.Weigh.g.Time3 - data$Dry.Weigh.g.Time2)/data$Dry.Weigh.g.Time2)*100)/data$Days.Time3 #calculate growth rate for time2
# data$time4.growth <- (((data$Dry.Weigh.g.Time4 - data$Dry.Weigh.g.Time3)/data$Dry.Weigh.g.Time3)*100)/data$Days.Time4 #calculate growth rate for time2
# data$time5.growth <- (((data$Dry.Weigh.g.Time5 - data$Dry.Weigh.g.Time4)/data$Dry.Weigh.g.Time4)*100)/data$Days.Time5 #calculate growth rate for time2
# data$time6.growth <- (((data$Dry.Weigh.g.Time6 - data$Dry.Weigh.g.Time5)/data$Dry.Weigh.g.Time5)*100)/data$Days.Time6 #calculate growth rate for time2
# 
# range(na.omit(data$time1.growth))
# range(na.omit(data$time2.growth))
# range(na.omit(data$time3.growth))
# range(na.omit(data$time4.growth))
# range(na.omit(data$time5.growth))
# range(na.omit(data$time6.growth))
# 
# mn <- -0.09
# mx <- 0.33
# 
# data$time1.growth[data$time1.growth < mn] <- NA
# data$time2.growth[data$time2.growth < mn] <- NA
# data$time3.growth[data$time3.growth < mn] <- NA
# data$time4.growth[data$time4.growth < mn] <- NA
# data$time5.growth[data$time5.growth < mn] <- NA
# data$time6.growth[data$time6.growth < mn] <- NA
# 
# data$time1.growth[data$time1.growth > mx] <- NA
# data$time2.growth[data$time2.growth > mx] <- NA
# data$time3.growth[data$time3.growth > mx] <- NA
# data$time4.growth[data$time4.growth > mx] <- NA
# data$time5.growth[data$time5.growth > mx] <- NA
# data$time6.growth[data$time6.growth > mx] <- NA
# 
# #reduced dataset of mass normalized growth
# G.mass <- data[,c(1:5,34:39)]
# 
# mean1 <- aggregate(time1.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se1 <- aggregate(time1.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T1 <- cbind(mean1, se1$time1.growth)
# colnames(T1) <- c("Parental.Performance", "Treatment",	"mean", "se")
# T1$Time <- "April 01"
# 
# mean2 <- aggregate(time2.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se2 <- aggregate(time2.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T2 <- cbind(mean2, se2$time2.growth)
# colnames(T2) <- c("Parental.Performance", "Treatment",  "mean", "se")
# T2$Time <- "April 20"
# 
# mean3 <- aggregate(time3.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se3 <- aggregate(time3.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T3 <- cbind(mean3, se3$time3.growth)
# colnames(T3) <- c("Parental.Performance", "Treatment",  "mean", "se")
# T3$Time <- "May 25"
# 
# mean4 <- aggregate(time4.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se4 <- aggregate(time4.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T4 <- cbind(mean4, se4$time4.growth)
# colnames(T4) <- c("Parental.Performance", "Treatment",  "mean", "se")
# T4$Time <- "June 17"
# 
# mean5 <- aggregate(time5.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se5 <- aggregate(time5.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T5 <- cbind(mean5, se5$time5.growth)
# colnames(T5) <- c("Parental.Performance", "Treatment",  "mean", "se")
# T5$Time <- "July 15"
# 
# mean6 <- aggregate(time6.growth ~ Parental.Performance*Treatment, data=G.mass, mean, na.rm=T)
# se6 <- aggregate(time6.growth ~ Parental.Performance*Treatment, data=G.mass, std.error, na.rm=T)
# T6 <- cbind(mean6, se6$time6.growth)
# colnames(T6) <- c("Parental.Performance", "Treatment",  "mean", "se")
# T6$Time <- "August 12"
# 
# 
# G <- rbind(T1, T2, T3, T4, T5, T6)
# G$TS <- c("BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH")
# 
# 
# 
# FigX <- ggplot(G, aes(x=Time, y=mean, group=TS)) + 
#   geom_errorbar(aes(ymin=G$mean-G$se, ymax=G$mean+G$se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
#   geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
#   scale_shape_manual(values=c(1,19)) + #sets point shape manually
#   geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
#   xlab("Time") + #Label the X Axis
#   ylab("Growth % per Day") + #Label the Y Axis
#   #ylim(-0.05, 0.3) + #set Y limits
#   theme_bw() + #Set the background color
#   theme(axis.line = element_line(color = 'black'), #Set the axes color
#         axis.title=element_text(size=14,face="bold"), #Set axis format
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background =element_blank(), #Set the plot background
#         legend.key = element_blank()) #Set plot legend key
# FigX

#####Growth Normalized to Surface area#####
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

data.sa <- reshape(BW.data, timevar = "TimePoint", drop = c("Tank", "QC", "Mass.g","Temp.C"), idvar=c("Fragment.ID","Parental.Performance", "Treatment", "Parent.ID","Cohort" ), direction="wide") #shape into wide format
data.sa <- merge(data.sa, SA.Data, by="Fragment.ID") #merge with the surface area file

#Growth Normalized to surface area mg CaCO3 cm-2 d-1
data.sa$time1.growth <- (((data.sa$Dry.Weigh.g.Time1 - data.sa$Dry.Weigh.g.Time0)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time1 #calculate growth rate for time1
data.sa$time2.growth <- (((data.sa$Dry.Weigh.g.Time2 - data.sa$Dry.Weigh.g.Time1)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time2 #calculate growth rate for time2
data.sa$time3.growth <- (((data.sa$Dry.Weigh.g.Time3 - data.sa$Dry.Weigh.g.Time2)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time3 #calculate growth rate for time3
data.sa$time4.growth <- (((data.sa$Dry.Weigh.g.Time4 - data.sa$Dry.Weigh.g.Time3)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time4 #calculate growth rate for time4
data.sa$time5.growth <- (((data.sa$Dry.Weigh.g.Time5 - data.sa$Dry.Weigh.g.Time4)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time5 #calculate growth rate for time5
data.sa$time6.growth <- (((data.sa$Dry.Weigh.g.Time6 - data.sa$Dry.Weigh.g.Time5)*1000/data.sa$SA.total.initial.cm2))/data.sa$Days.Time6 #calculate growth rate for time6

#view ranges to look for outliers 
range(na.omit(data.sa$time1.growth))
range(na.omit(data.sa$time2.growth))
range(na.omit(data.sa$time3.growth))
range(na.omit(data.sa$time4.growth))
range(na.omit(data.sa$time5.growth))
range(na.omit(data.sa$time6.growth))

#view histograms of each time points growth data
par(mfrow = c(3, 2)) # 3 row x 2 columns on one plot
hist(data.sa$time1.growth)
hist(data.sa$time2.growth)
hist(data.sa$time3.growth)
hist(data.sa$time4.growth)
hist(data.sa$time5.growth)
hist(data.sa$time6.growth)

#set low and upper bounds from ranges
mn <- -0.62
mx <- 1.2

#replace values outside min and max with NA
data.sa$time1.growth[data.sa$time1.growth < mn] <- NA
data.sa$time2.growth[data.sa$time2.growth < mn] <- NA
data.sa$time3.growth[data.sa$time3.growth < mn] <- NA
data.sa$time4.growth[data.sa$time4.growth < mn] <- NA
data.sa$time5.growth[data.sa$time5.growth < mn] <- NA
data.sa$time6.growth[data.sa$time6.growth < mn] <- NA

data.sa$time1.growth[data.sa$time1.growth > mx] <- NA
data.sa$time2.growth[data.sa$time2.growth > mx] <- NA
data.sa$time3.growth[data.sa$time3.growth > mx] <- NA
data.sa$time4.growth[data.sa$time4.growth > mx] <- NA
data.sa$time5.growth[data.sa$time5.growth > mx] <- NA
data.sa$time6.growth[data.sa$time6.growth > mx] <- NA

par(mfrow = c(3, 2)) # 3 row x 2 columns on one plot
hist(data.sa$time1.growth)
hist(data.sa$time2.growth)
hist(data.sa$time3.growth)
hist(data.sa$time4.growth)
hist(data.sa$time5.growth)
hist(data.sa$time6.growth)

#reduced dataset of surface area normalized growth
G.SA <- data.sa[,c(1:5,39:44)]

mean1 <- aggregate(time1.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se1 <- aggregate(time1.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T1 <- cbind(mean1, se1$time1.growth)
colnames(T1) <- c("Parental.Performance", "Treatment",  "mean", "se")
T1$Time <- "Time1"

mean2 <- aggregate(time2.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se2 <- aggregate(time2.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T2 <- cbind(mean2, se2$time2.growth)
colnames(T2) <- c("Parental.Performance", "Treatment",  "mean", "se")
T2$Time <- "Time2"

mean3 <- aggregate(time3.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se3 <- aggregate(time3.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T3 <- cbind(mean3, se3$time3.growth)
colnames(T3) <- c("Parental.Performance", "Treatment",  "mean", "se")
T3$Time <- "Time3"

mean4 <- aggregate(time4.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se4 <- aggregate(time4.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T4 <- cbind(mean4, se4$time4.growth)
colnames(T4) <- c("Parental.Performance", "Treatment",  "mean", "se")
T4$Time <- "Time4"

mean5 <- aggregate(time5.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se5 <- aggregate(time5.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T5 <- cbind(mean5, se5$time5.growth)
colnames(T5) <- c("Parental.Performance", "Treatment",  "mean", "se")
T5$Time <- "Time5"

mean6 <- aggregate(time6.growth ~ Parental.Performance*Treatment, data=G.SA, mean, na.rm=T)
se6 <- aggregate(time6.growth ~ Parental.Performance*Treatment, data=G.SA, std.error, na.rm=T)
T6 <- cbind(mean6, se6$time6.growth)
colnames(T6) <- c("Parental.Performance", "Treatment",  "mean", "se")
T6$Time <- "Time6"

Growth.SA <- rbind(T1, T2, T3, T4, T5, T6)
Growth.SA$TS <- c("BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH","BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH", "BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH", "BL-AMB", "UNBL-AMB", "BL-HIGH", "UNBL-HIGH")

Fig.G.SA <- ggplot(Growth.SA, aes(x=Time, y=mean, group=TS)) + 
  geom_errorbar(aes(ymin=Growth.SA$mean-Growth.SA$se, ymax=Growth.SA$mean+Growth.SA$se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #plot sem
  geom_point(aes(shape=Parental.Performance), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(1,19)) + #sets point shape manually
  geom_line(aes(linetype=Treatment), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #Label the X Axis
  ylab("Growth mg cm-2 d-1") + #Label the Y Axis
  scale_x_discrete(labels=c("April01","April20","May25","June17", "July15", "August12")) +
  theme_bw() + #Set the background color
  theme(axis.line = element_line(color = 'black'), #Set the axes color
        axis.title=element_text(size=14,face="bold"), #Set axis format
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background =element_blank(), #Set the plot background
        legend.key = element_blank()) #Set plot legend key
Fig.G.SA

setwd("/Users/hputnam/MyProjects/Mcap_PGA_TGA/RAnalysis/Output/") #set working
ggsave(file="Calcification.pdf", Fig.G.SA, width = 11, height = 6, units = c("in"))



###Repeated Measures ANOVA
G.RM <- melt(G.SA) #reshape into long format

Growth.RM.lme <- lme(value ~ variable*Parental.Performance*Treatment, random = ~ variable|Fragment.ID, data=na.omit(G.RM)) #repeated measures ANOVA with random intercept but not slope (clonal fragments expect to respond the same)
Growth.results <- summary(Growth.RM.lme) #view RM ANOVA summary
Growth.stats <- anova(Growth.RM.lme) #view F and p values
Growth.stats

Growth.RM.posthoc <- lsmeans(Growth.RM.lme, specs=c("variable","Parental.Performance")) #calculate MS means
Growth.RM.posthoc #view results
Growth.RM.posthoc.p <- contrast(Growth.RM.posthoc, method="pairwise", by=c("variable")) #contrast treatment groups within a species at each time point
Growth.RM.posthoc.p #view results

###Testing ANOVA Assumptions
###Data are normally distributed
hist(Growth.RM.lme$residuals) #histogram
qqnorm(Growth.RM.lme$residuals) # normal quantile plot

