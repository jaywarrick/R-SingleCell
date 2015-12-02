# rm(list = ls())

# xmin=-100
# xmax=500000
# ymin=-100
# ymax=500000
# xThresh=74.79895
# yThresh=168.0746
# xCross=0
# yCross=0.002
# xTransition=xThresh
# yTransition=yThresh
# xLinLogRatio=200
# yLinLogRatio=xLinLogRatio
# xTicks=c(0,100,1000,10000,100000,1000000)
# yTicks=xTicks

source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/logicle.R')
source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/plotFACS.R')
source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/getSingleStats.R')
source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/getDoubleStats.R')
source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/drawStats.R')
source('/Users/jaywarrick/Google Drive/SingleCell/FACS Data/R FACS Plotting Functions/drawLogicleAxis.R')
source('/Users/jaywarrick/Google Drive/SingleCell/Figures/plotHelperFunctions.R')
library(flowCore)
library(flowViz)
library(MASS)
library(RColorBrewer)
library(foreign)
library(hash)

microscopeFolder <- '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/'

# Move to the folder where Microwell plots are kept
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/Microwell/')

results <- list()

data1 <- read.table(paste0(microscopeFolder, 'MasterDataFile.txt'), header=TRUE)
WT0 <- data.frame(G=subset(data1, Time=='0')$G.Final, R=subset(data1, Time=='0')$R.Final)
medWT <- c(median(WT0$G), median(WT0$R))
madWT <- c(mad(WT0$G), mad(WT0$R))
threshWT <- medWT + 3*madWT


times <- unique(data1$Time)
results1 <- list()
for(t in times)
{
    temp <- data.frame(G=subset(data1, Time==t & Measurement=='G')$Value, R=subset(data1, Time==t & Measurement=='R')$Value)
    stats <- getDoubleStats(temp$G, temp$R, xThresh=threshWT[1], xCross=0, yThresh=threshWT[2], yCross=0.004)
    results1[[t]] <- stats
}


data2 <- read.arff(paste(microscopeFolder, 'x1_y1.arff', sep=''))
M51R0 <- data.frame(G=subset(data2, Time=='0' & Measurement=='G')$Value, R=subset(data2, Time=='0' & Measurement=='R')$Value)
medM51R <- c(median(M51R0$G), median(M51R0$R))
madM51R <- c(mad(M51R0$G), mad(M51R0$R))
threshM51R <- medM51R + 3*madM51R


times <- unique(data2$Time)
results2 <- list()
for(t in times)
{
    temp <- data.frame(G=subset(data2, Time==t & Measurement=='G')$Value, R=subset(data2, Time==t & Measurement=='R')$Value)
    stats <- getDoubleStats(temp$G, temp$R, xThresh=threshM51R[1], xCross=0, yThresh=threshM51R[2], yCross=0.004)
    results2[[t]] <- stats
}

data3 <- read.arff(paste(microscopeFolder, 'x1_y2.arff', sep=''))
MOCK0 <- data.frame(G=subset(data2, Time=='0' & Measurement=='G')$Value, R=subset(data2, Time=='0' & Measurement=='R')$Value)
medMOCK <- c(median(MOCK0$G), median(MOCK0$R))
madMOCK <- c(mad(MOCK0$G), mad(MOCK0$R))
threshMOCK <- medMOCK + 3*madMOCK
print(threshMOCK)


times <- unique(data3$Time)
results3 <- list()
for(t in times)
{
    temp <- data.frame(G=subset(data3, Time==t & Measurement=='G')$Value, R=subset(data3, Time==t & Measurement=='R')$Value)
    stats <- getDoubleStats(temp$G, temp$R, xThresh=threshMOCK[1], xCross=0, yThresh=threshMOCK[2], yCross=0.004)
    results3[[t]] <- stats
}


######## Plot 1 ###############
fun <- function(hashObj, hashKey)
{
    return(hashObj[[hashKey]])
}
times2 <- (10/60)*as.numeric(as.character(times))+2.833333
pdf('WT_Percentages_Plot.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(min(times2),max(times2)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
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
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

################ Plot 2 ###############
times2 <- (10/60)*as.numeric(as.character(times))+2.833333
pdf('M51R_Percentages_Plot.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(min(times2),max(times2)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
lines(times2, lapply(results2, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times2, lapply(results2, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(results2, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(results2, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(results2, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex = 0.75)
points(times2, lapply(results2, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex = 0.75)
points(times2, lapply(results2, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex = 0.75)
points(times2, lapply(results2, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex = 0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

################ Plot 3 ###############
times3 <- (10/60)*as.numeric(as.character(times))+2.833333
pdf('MOCK_Percentages_Plot.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(min(times3),max(times3)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
lines(times3, lapply(results3, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times3, lapply(results3, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times3, lapply(results3, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times3, lapply(results3, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times3, lapply(results3, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex = 0.75)
points(times3, lapply(results3, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex = 0.75)
points(times3, lapply(results3, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex = 0.75)
points(times3, lapply(results3, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex = 0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()


