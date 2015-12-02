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

source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/logicle.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/plotFACS.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/getSingleStats.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/getDoubleStats.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/drawStats.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/drawLogicleAxis.R')
source('/Users/jaywarrick/Public/DropBox/Github/R-SingleCell/plotHelperFunctions.R')
library(flowCore)
library(flowViz)
library(MASS)
library(RColorBrewer)
library(foreign)
library(hash)

FACSFolder <- '/Users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data/FCS Data/Compiled/'
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figures/Plots/FC')

# Crossover for FACS = 0.002
linLogRatio <- 300
ticks <- c(10,100,1000,10000,100000,1000000)
tickLabels <- expression(10^1,10^2,10^3,10^4,10^5,10^6)
plotFACSParams <- list(title='title', temp=c(), xmin=-100, xmax=200000, ymin=-200, ymax=300000, xCross=0, yCross=0.002, xLinLogRatio=linLogRatio, yLinLogRatio=linLogRatio, xTicks=ticks, yTicks=ticks, xTickLabels=tickLabels, yTickLabels=tickLabels, time='-1')


thresh <- list()
med <- list()
files <- c('WT_6_B.fcs', 'WT_12_C.fcs', 'WT_S_A.fcs', 'M51R_6_A.fcs', 'M51R_12_C.fcs', 'M51R_S_B.fcs', 'Mock_6_C.fcs', 'Mock_12_C.fcs', 'Mock_S_B.fcs')
newNames <- c('FC_WT_6','FC_WT_12','FC_WT_18','FC_M51R_6','FC_M51R_12','FC_M51R_18','FC_MOCK_6','FC_MOCK_12','FC_MOCK_18')
# files <- c('Mock_6_C.fcs', 'Mock_12_C.fcs', 'Mock_S_B.fcs')

data <- read.FCS(paste(FACSFolder, 'Mock_6_C.fcs', sep=''), transformation=FALSE)
temp <- data.frame(G=exprs(data$'488 B 530/30-A')[,1], R=exprs(data$'561 D 610/20-A')[,1])
med$Early <- c(median(temp$G), median(temp$R))
mad <- c(mad(temp$G), mad(temp$R))
thresh$Early <- 3*mad

data <- read.FCS(paste(FACSFolder, 'Mock_12_C.fcs', sep=''), transformation=FALSE)
temp <- data.frame(G=exprs(data$'488 B 530/30-A')[,1], R=exprs(data$'561 D 610/20-A')[,1])
med$Middle <- c(median(temp$G), median(temp$R))
mad <- c(mad(temp$G), mad(temp$R))
thresh$Middle <- 3*mad

data <- read.FCS(paste(FACSFolder, 'Mock_S_B.fcs', sep=''), transformation=FALSE)
temp <- data.frame(G=exprs(data$'488 B 530/30-A')[,1], R=exprs(data$'561 D 610/20-A')[,1])
med$Late <- c(median(temp$G), median(temp$R))
mad <- c(mad(temp$G), mad(temp$R))
thresh$Late <- 3*mad

for(f in files)
{
    data <- read.FCS(paste(FACSFolder, f, sep=''), transformation=FALSE)
    temp <- data.frame(G=exprs(data$'488 B 530/30-A')[,1], R=exprs(data$'561 D 610/20-A')[,1])
    if(length(grep('6', f, fixed=TRUE)) > 0)
    {
        gateX <- thresh$Early[1]
        gateY <- thresh$Early[2]
        temp$G <- temp$G - med$Early[1]
        temp$R <- temp$R - med$Early[2]
    }
    if(length(grep('12', f, fixed=TRUE)) > 0)
    {
        gateX <- thresh$Middle[1]
        gateY <- thresh$Middle[2]
        temp$G <- temp$G - med$Middle[1]
        temp$R <- temp$R - med$Middle[2]
    }
    if(length(grep('18', f, fixed=TRUE)) > 0)
    {
        gateX <- thresh$Late[1]
        gateY <- thresh$Late[2]
        temp$G <- temp$G - med$Late[1]
        temp$R <- temp$R - med$Late[2]
    }
    if(length(grep('S', f, fixed=TRUE)) > 0)
    {
        gateX <- thresh$Late[1]
        gateY <- thresh$Late[2]
        temp$G <- temp$G - med$Late[1]
        temp$R <- temp$R - med$Late[2]
    }
    print(gateX)
    print(gateY)
    print(f)
    pdf(paste(newNames[which(files==f)],'.pdf',sep=''),width=1.1*W,height=1.1*H)
    plotFACSParams <- merge.lists(plotFACSParams, list(title=f, temp=temp, xThresh=gateX, yThresh=gateY, xTransition=gateX, yTransition=gateY))
    do.call(plotFACS, plotFACSParams)
    dev.off()
}

# Plot the summary plots
times2 <- c(6,12,18)
results1 <- list()
for(f in files)
{
     data <- read.FCS(paste(FACSFolder, f, sep=''), transformation=FALSE)
     temp <- data.frame(G=exprs(data$'488 B 530/30-A')[,1], R=exprs(data$'561 D 610/20-A')[,1])
     if(length(grep('6', f, fixed=TRUE) > 0))
     {
          gateX <- thresh$Early[1]
          gateY <- thresh$Early[2]
          temp$G <- temp$G - med$Early[1]
          temp$R <- temp$R - med$Early[2]
     }
     if(length(grep('12', f, fixed=TRUE) > 0))
     {
          gateX <- thresh$Middle[1]
          gateY <- thresh$Middle[2]
          temp$G <- temp$G - med$Middle[1]
          temp$R <- temp$R - med$Middle[2]
     }
     if(length(grep('S', f, fixed=TRUE) > 0))
     {
          gateX <- thresh$Late[1]
          gateY <- thresh$Late[2]
          temp$G <- temp$G - med$Late[1]
          temp$R <- temp$R - med$Late[2]
     }
     stats <- getDoubleStats(temp$G, temp$R, xThresh=gateX, xCross=0, yThresh=gateY, yCross=0.002)
     results1[[f]] <- stats
}


######## Plot 1 ###############
fun <- function(hashObj, hashKey)
{
     return(hashObj[[hashKey]])
}
pdf('FC_WT_Summary.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(2.8333,max(times2)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
resultsWT <- results1[1:3]
lines(times2, lapply(resultsWT, fun, hashKey='XY % ++'), type='l', col='goldenrod3', lwd=2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(resultsWT, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex=0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

######## Plot M51R ###############
fun <- function(hashObj, hashKey)
{
     return(hashObj[[hashKey]])
}
paperParams(1,1)
pdf('FC_M51R_Summary.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(2.8333,max(times2)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
resultsWT <- results1[4:6]
lines(times2, lapply(resultsWT, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(resultsWT, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex=0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()

######## Plot MOCK ###############
fun <- function(hashObj, hashKey)
{
     return(hashObj[[hashKey]])
}
paperParams(1,1)
pdf('FC_MOCK_Summary.pdf',width=1.1*W,height=1.1*H)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(xlimit=c(2.8333,max(times2)), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', plotaxes=FALSE, xcol='black', ycol='black')
resultsWT <- results1[7:9]
lines(times2, lapply(resultsWT, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2)
lines(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2)
points(times2, lapply(resultsWT, fun, hashKey='XY % ++'), col='goldenrod3', lwd=2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % +-'), col = 'darkgreen', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % -+'), col = 'darkred', lwd = 2, cex=0.75)
points(times2, lapply(resultsWT, fun, hashKey='XY % --'), col = 'black', lwd = 2, cex=0.75)
axis(1, at=c(3,6,12,18))
axis(2)
box(col='black', lwd=1.5)
legend('left',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
dev.off()


