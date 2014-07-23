rm(list=ls())
library(drc)
library(zoo)
library(plyr)
library(data.table)
source('/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/postprocessingHelperFunctions.R')

globalFctDefA <- list(name='global',fct=fctA, ssfct=ssfctA, names=c('alpha','tau'));

cont <- drmc(constr = FALSE, errorm = TRUE, maxIt = 1000, method="SANN", 
             noMessage = FALSE, relTol = 1e-07, rmNA=FALSE, useD = FALSE, 
             trace = TRUE, otrace = TRUE, warnVal = -1, dscaleThres = 1e-15, rscaleThres = 1e-15)

# Read the m51R data
data <- read.table('/users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/M51R_Data_1Cell.txt', header=TRUE)
IDs <- unique(data$ID)

# For each microwell, if we have both green and red signal, and R^2 fits of > 0.9, and alpha's greater than 0,
# store the alpha for Red and Green along with detections times and max levels
#
# Create a table for storing calculated information
# params <- data.frame(Cell=numeric(0), Max.R=numeric(0), Max.G=numeric(0), Time.Delay.R=numeric(0), Time.Delay.G=numeric(0), Time.Rise.R=numeric(0), Time.Rise.G=numeric(0), Time.Peak.R=numeric(0), Time.Peak.G=numeric(0), Time.Death=numeric(0), Mod0.Alpha.R=numeric(0), Mod0.Alpha.G=numeric(0), Mod0.R2.R=numeric(0), Mod0.R2.G=numeric(0), Mod0.Flag.Quality=numeric(0)))
params <- data.frame();
newDataR <- data.frame();
newDataG <- data.frame();
results <- list();
for(id in IDs)
{
    temp <- subset(data, ID==id)
    Int.Max.R <- max(temp$R)
    Int.Max.G <- max(temp$G)
    if(Int.Max.R > 0 & Int.Max.G > 0) # RFP+ GFP+
    {
        results <- getDualParams(id, temp, Int.Max.R, Int.Max.G)
        newDataR <- rbind(newDataR, results$newDataR)
        newDataG <- rbind(newDataG, results$newDataG)
        params <- rbind(params, results$newParams)
    } else if(Int.Max.R > 0) # RFP+ GFP-
    {
        results <- getSingleParams(id, temp, Int.Max.R)
        newDataR <- rbind(newDataR, results$newDataR)
        params <- rbind(params, results$newParams)
    }
}
params$Time.Delta.Delay <- params$Time.Delay.G-params$Time.Delay.R
params$Time.Rise.R <- params$Time.Peak.R-params$Time.Delay.R
params$Time.Rise.G <- params$Time.Peak.G-params$Time.Delay.G
deaths <- !is.na(params$Time.Death)
params$Flag.Poor.Death.Detection <- NA
params$Flag.Poor.Death.Detection[deaths] <- (params$Time.Death[deaths] < params$Time.Peak.R[deaths]) | (params$Time.Death[deaths] < params$Time.Max.R[deaths]) | (params$Time.Death[deaths] < params$Time.Peak.G[deaths]) | (params$Time.Death[deaths] < params$Time.Max.G[deaths])
params$NLFit.Alpha.R <- NA
params$NLFit.Alpha.G <- NA
params$NLFit.Tau.R <- NA
params$NLFit.Tau.G <- NA
params$NLFit.R2.R <- NA
params$NLFit.R2.G <- NA
params <- data.table(params)
params[,Flag.R:=data[data$ID==ID,]$Flag.R[1], by=ID]
params[,Flag.G:=data[data$ID==ID,]$Flag.G[1], by=ID]


# Plot the 4 fitted timepoints for all data to get a sense check
plot(c(),c(), xlim=c(min(newDataR$time),max(newDataR$time)), ylim=c(1,3000), log='y')
duh <- data.table(newDataR)
duh[,lines(time,R), by=ID]

# Calculate the weights for the fitting function (use Poisson model of noise in microscopy data)
BG <- 350
newDataR$w <- 1/((0.059^2)*(newDataR$R+BG))
newDataG$w <- 1/((0.059^2)*(newDataG$G+BG))


# Fit the NON-LOGGED data with an EXPONENTIAL fit with CONSTANT delay
subsetDataR <- subset(newDataR, ID %in% subset(params, LMFit.Flag.Quality.R==1)$ID & Flag.R==TRUE) # (RFP+ GFP+) & (RFP+ GFP-)
subsetDataG <- subset(newDataG, ID %in% subset(params, LMFit.Flag.Quality.G==1)$ID & Flag.R==TRUE & Flag.G==TRUE) # (RFP+ GFP+) only

# Fitting non-logged data. A means 1 single tau for all curves, 2 means with variable tau, R means red and G means green
print('fitting modA.R')
modA.R <- drm(R~time,data=subsetDataR,ID,pmodels=data.frame(ID,1),fct=globalFctDefA, weights=subsetDataR$w, na.action=na.omit, separate=FALSE, lowerl=c(0,0), upperl=c(4,50))
summary(modA.R)
plot(modA.R, ylim=c(1,1000), log='y')
print('fitting modA.G')
modA.G <- drm(G~time,data=subsetDataG,ID,pmodels=data.frame(ID,1),fct=globalFctDefA, weights=subsetDataG$w, na.action=na.omit, separate=FALSE, lowerl=c(0,0))
summary(modA.G)
plot(modA.G, ylim=c(1,1000), log='y')


# VARIABLE delay (tau)
print('fitting modB.R')
modB.R <- drm(R~time,data=subsetDataR,ID,pmodels=data.frame(ID,ID),fct=globalFctDefA, weights=subsetDataR$w, na.action=na.omit, separate=FALSE, lowerl=c(0,0))
summary(modB.R)
subsetDataR$fitted <- fitted(modB.R)
subsetDataR$residuals <- residuals(modB.R)
plot(modB.R, ylim=c(1,1000), log='y')
print('fitting modB.G')
modB.G <- drm(G~time,data=subsetDataG,ID,pmodels=data.frame(ID,ID),fct=globalFctDefA, weights=subsetDataG$w, na.action=na.omit, separate=FALSE, lowerl=c(0,0))
summary(modB.G)
subsetDataG$fitted <- fitted(modB.G)
subsetDataG$residuals <- residuals(modB.G)
plot(modB.G, ylim=c(1,1000), log='y')

params$NLFit.Flag.Quality.R <- NA
params$NLFit.Flag.Quality.G <- NA
for(id in unique(subsetDataR$ID)) # (RFP+ GFP+) & (RFP+ GFP-)
{
    print(id)
    temp <- subset(data, ID==id)
    pID <- which(params$ID==id)
    
    # get red NLFit information
    R <- NA
    R <- subset(subsetDataR, ID==id)
    params$NLFit.Alpha.R[pID] <- getPar(id, 'alpha', modB.R)
    params$NLFit.Tau.R[pID] <- getPar(id, 'tau', modB.R)
    R2.R <- getWeightedR2_Vectors(y=R$R, r=R$residuals, f=R$fitted, w=R$w)
    params$NLFit.R2.R[pID] <- R2.R
    params$NLFit.Flag.Quality.R[pID] <- (params$NLFit.Alpha.R[pID] > 0) & (params$NLFit.R2.R[pID] >= 0.9)
    
    # get green NLFit information if relevant (i.e. RFP+ GFP+ only)
    G <- NA;
    if(params[pID]$LMFit.Flag.Quality.G==TRUE & (params[pID]$Flag.G > 0))
    {
        hasGreen <- TRUE;
        G <- subset(subsetDataG, ID==id)
        params$NLFit.Alpha.G[pID] <- getPar(id, 'alpha', modB.G)
        params$NLFit.Tau.G[pID] <- getPar(id, 'tau', modB.G)
        R2.G <- getWeightedR2_Vectors(y=G$G, r=G$residuals, f=G$fitted, w=G$w)
        params$NLFit.R2.G[pID] <- R2.G
        params$NLFit.Flag.Quality.G[pID] <- (params$NLFit.Alpha.G[pID] > 0) & (params$NLFit.R2.G[pID] >= 0.9)
    }
    
    # plot data and save the NLFit information to params
    plotID(params[pID], temp, R, G);
}
params$ID[which(params$NLFit.Flag.Quality.R==FALSE)]
params$ID[which(params$NLFit.Flag.Quality.G==FALSE)]

anova(modA.R,modB.R)
anova(modA.G,modB.G)

library(RCulr)
source("http://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt")
source_https("https://raw.github.com/talgalili/R-code-snippets/master/siegel.tukey.r")
duh <- siegel.tukey(params$Time.Delay.R, params$Time.Delta.Delay[!is.na(params$Time.Delta.Delay)], adjust.median=TRUE)

write.table(params, file='/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/TrajectoryParameters.txt')
params2 <- read.table('/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/TrajectoryParameters.txt', header=TRUE)

# Include assessment of Death (0 = No Plateau Verified, 1 = Plateau, 2 = Lysis, 3 = Neither)
deathTable <- read.table('/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/Files from Andrea/M51RRed_10.19.13.txt', sep='\t', header=TRUE)
deathTable <- rbind(deathTable, read.table('/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/Files from Andrea/M51RBoth_10.19.13.txt', sep='\t', header=TRUE))
deathTable$Image2 <- deathTable$Image + 36
deathTable$ID <- deathTable$Image2*1000 + deathTable$ROI
params2 <- merge(params2, deathTable[,c('ID','Death')])
write.table(params2, file='/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/TrajectoryParameters2.txt')

write.table(data, file='/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/MasterTable_M51R.txt')

theTimes <- unique(data$time)
theTime <- theTimes[which(theTimes > 17.5)[1]]

fplotData <- subset(data, time==theTime)
plot(fplotData$G,fplotData$R, log='xy')