# File for plotting the PDC plots for microwells for M51R and N1
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/logicle.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/plotFACS.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/getSingleStats.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/getDoubleStats.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/drawStats.R')
source('/Users/jaywarrick/Public/DropBox/GitHub/R-SingleCell/drawLogicleAxis.R')
source('/Users/jaywarrick/Public/Dropbox/GitHub/R-SingleCell/plotHelperFunctions.R')
library(flowCore)
library(flowViz)
library(MASS)
library(RColorBrewer)
library(foreign)
library(hash)
library(traitr)

# Set plotting parameters
# Crossover for Yin Scope = 0.0007
linLogRatio <- 35
ticks <- c(0,10,100,1000,10000,100000,1000000)
tickLabels <- expression(0,10^1,10^2,10^3,10^4,10^5,10^6)
plotFACSParams <- list(title='title', temp=c(), xmin=-5, xmax=4000, ymin=-30, ymax=11000, xCross=0, yCross=0.00007, xLinLogRatio=linLogRatio, yLinLogRatio=linLogRatio, xTicks=ticks, yTicks=ticks, xTickLabels=tickLabels, yTickLabels=tickLabels, xcol='darkgreen', ycol='darkred', time='-1')

########### M51R DC Plots ##############

# Store the timepoints that we'll creat static plots
results <- list()
data <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Data_1Cell.txt', header=TRUE)
data <- subset(data, Flag.R)
results$M0 <- data.frame(G=subset(data, Time=='0')$G.Final, R=subset(data, Time=='0')$R.Final)
results$M1 <- data.frame(G=subset(data, Time=='9')$G.Final, R=subset(data, Time=='9')$R.Final)
results$M2 <- data.frame(G=subset(data, Time=='21')$G.Final, R=subset(data, Time=='21')$R.Final)
results$M3 <- data.frame(G=subset(data, Time=='45')$G.Final, R=subset(data, Time=='33')$R.Final)
medM <- c(median(results$M0$G), median(results$M0$R))
madM <- c(mad(results$M0$G), mad(results$M0$R))
threshM <- medM + 3*madM
print(threshM)

# Make the M51R movie files
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell/DC_M51R')
for(t in unique(data$Time))
{
     temp <- data.frame(G=subset(data, Time==t)$G.Final, R=subset(data, Time==t)$R.Final)
     pdf(paste('Time_',sprintf('%03.f',as.numeric(t)),'.pdf',sep=''),width=1.1*W,height=1.1*H)
     paperParams(1,1)
     tempT <- data[data$Time==t,'time'][1]
     dTime <- sprintf('%0.2f',tempT)
     plotFACSParams <- merge.lists(plotFACSParams, list(title='dummy', temp=temp, xThresh=threshM[1], yThresh=threshM[2], xTransition=threshM[1], yTransition=threshM[2], time=dTime))
     do.call(plotFACS, plotFACSParams)
     dev.off()
}
cmd1 <- 'convert -verbose -background white -density 300 *.pdf -quality 100 %03d.jpg'
cmd2 <- 'ffmpeg -y -r 6 -i %03d.jpg -vcodec mjpeg -qscale 1 -an output.avi'
system(cmd1)
system(cmd2)

# Make the M51R summary plots
times <- unique(data$Time)
results1 <- list()
for(t in times)
{
     temp <- data.frame(G=subset(data, Time==t)$G.Final, R=subset(data, Time==t)$R.Final)
     stats <- getDoubleStats(temp$G, temp$R, xThresh=threshM[1], xCross=0, yThresh=threshM[2], yCross=0.004)
     results1[[as.character(t)]] <- stats
}
fun <- function(hashObj, hashKey)
{
     return(hashObj[[hashKey]])
}
times2 <- unique(data$time)
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell')
pdf('Microwell_M51R_Summary.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(3,18), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
lines(times2, lapply(results1, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times2, lapply(results1, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(results1, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(results1, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(results1, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex = 0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('top',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

########### WT (N1) DC Plots ##############

# Store the timepoints that we'll creat static plots
data <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Data_1Cell.txt', header=TRUE)
data <- subset(data, Flag.R)
results$N0 <- data.frame(G=subset(data, Time=='0')$G.Final, R=subset(data, Time=='0')$R.Final)
results$N1 <- data.frame(G=subset(data, Time=='9')$G.Final, R=subset(data, Time=='9')$R.Final)
results$N2 <- data.frame(G=subset(data, Time=='21')$G.Final, R=subset(data, Time=='21')$R.Final)
results$N3 <- data.frame(G=subset(data, Time=='45')$G.Final, R=subset(data, Time=='33')$R.Final)
medN <- c(median(results$N0$G), median(results$N0$R))
madN <- c(mad(results$N0$G), mad(results$N0$R))
threshN <- medN + 3*madN
print(threshN)

# Make the WT movie files
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell/DC_WT')
for(t in unique(data$Time))
{
     temp <- data.frame(G=subset(data, Time==t)$G.Final, R=subset(data, Time==t)$R.Final)
     pdf(paste('Time_',sprintf('%03.f',as.numeric(t)),'.pdf',sep=''),width=1.1*W,height=1.1*H)
     paperParams(1,1)
     tempT <- data[data$Time==t,'time'][1]
     dTime <- sprintf('%0.2f',tempT)
     plotFACSParams <- merge.lists(plotFACSParams, list(title='dummy', temp=temp, xThresh=threshN[1], yThresh=threshN[2], xTransition=threshN[1], yTransition=threshN[2], time=dTime))
     do.call(plotFACS, plotFACSParams)
     dev.off()
}
cmd1 <- 'convert -verbose -background white -density 300 *.pdf -quality 100 %03d.jpg'
cmd2 <- 'ffmpeg -y -r 6 -i %03d.jpg -vcodec mjpeg -qscale 1 -an output.avi'
system(cmd1)
system(cmd2)

# Make the WT summary plot
data <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Data_1Cell.txt', header=TRUE)
data <- subset(data, Flag.R)
times <- unique(data$Time)
results1 <- list()
for(t in times)
{
     temp <- data.frame(G=subset(data, Time==t)$G.Final, R=subset(data, Time==t)$R.Final)
     stats <- getDoubleStats(temp$G, temp$R, xThresh=threshM[1], xCross=0, yThresh=threshM[2], yCross=0.004)
     results1[[as.character(t)]] <- stats
}
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell')
pdf('Microwell_WT_Summary.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(3,18), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
lines(times2, lapply(results1, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times2, lapply(results1, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(results1, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(results1, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(results1, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex = 0.75)
points(times2, lapply(results1, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex = 0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('top',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

########### M51R and WT static plots ############

# Go to the directory for Microwell plots
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell')
items <- names(results)
for(item in items)
{
     if(length(grep('M',item,fixed=TRUE)) > 0)
     {
          # then it's M51R data
          pdf(paste('M51R_',item,'.pdf',sep=''),width=1.1*W,height=1.1*H)
          paperParams(1,1)
          plotFACSParams <- merge.lists(plotFACSParams, list(title='dummy', temp=results[[item]], xThresh=threshM[1], yThresh=threshM[2], xTransition=threshM[1], yTransition=threshM[2], time='-1'))
          do.call(plotFACS, plotFACSParams)
          dev.off()
     } else
     {
          # its' WT data
          pdf(paste('WT_',item,'.pdf',sep=''),width=1.1*W,height=1.1*H)
          paperParams(1,1)
          plotFACSParams <- merge.lists(plotFACSParams, list(title='dummy', temp=results[[item]], xThresh=threshN[1], yThresh=threshN[2], xTransition=threshN[1], yTransition=threshN[2], time='-1'))
          do.call(plotFACS, plotFACSParams)
          dev.off()
     }
}

