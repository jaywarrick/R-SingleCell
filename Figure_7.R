
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

# Have data object prepared (i.e. use preProcessing.R & postProcessing.R) or read in the appropriate microscope data
data <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Data_1Cell.txt', header=TRUE)
thresholds <- read.table('/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Thresholds.txt', header=TRUE)
setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Figure Masters/Plots/Microwell')

threshR <- thresholds$R
threshG <- thresholds$G

theTimes <- unique(data$Time) # 'Time' is frame number while 'time' is hpi
thetimes <- unique(data$time) # 'Time' is frame number while 'time' is hpi
keyFrames <- theTimes[c(which(thetimes > 6)[1]-1, which(thetimes > 12)[1]-1, which(thetimes > 18)[1]-1)] #  the frames near 6, 12, and 18 hrs
print(paste('keyFrames :', keyFrames))
fractions <- data.frame(t=numeric(0), pp=numeric(0), pm=numeric(0), mp=numeric(0), mm=numeric(0))
for(frame in unique(data$Time))
{
    temp <- subset(data, Time==frame & (Flag.R | Flag.G))
    pdf(paste('MicrowellFACS_N1_',temp$Time[1],'.pdf',sep=''),width=1.1*W,height=1.1*H)
    stats <- plotFACS(title='dummy', temp, xThresh=threshG, yThresh=threshR, xTransition=threshG, yTransition=threshR, xmax=4000, ymax=13000, yCross=0.0007)
    fractions <- rbind(fractions, data.frame(t=temp$time[1], pp=stats$'XY % ++', pm=stats$'XY % +-', mp=stats$'XY % -+', mm=stats$'XY % --'))
    dev.off()
    #     browser()
}

pdf('MicrowellSummary_N1.pdf', height=1.1*H, width=1.1*W)
paperParams(1,1)
par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(letter='',xlimit=c(3,18), ylimit=c(0,100), xlabel='Time [hpi]', ylabel='% of Population', xcol='black', ycol='black', xticks=c('3','6','12','18'))

lines(fractions$t, fractions$pp, pch=21, col='goldenrod3', cex=0.75, lwd=2)
lines(fractions$t, fractions$pm, pch=21, col='darkgreen', cex=0.75, lwd=2)
lines(fractions$t, fractions$mp, pch=21, col='darkred', cex=0.75, lwd=2)
lines(fractions$t, fractions$mm, pch=21, col='black', cex=0.75, lwd=2)
points(fractions$t, fractions$pp, pch=21, col='goldenrod3', cex=0.75, lwd=2)
points(fractions$t, fractions$pm, pch=21, col='darkgreen', cex=0.75, lwd=2)
points(fractions$t, fractions$mp, pch=21, col='darkred', cex=0.75, lwd=2)
points(fractions$t, fractions$mm, pch=21, col='black', cex=0.75, lwd=2)
legend('top',c('[+,+]','[+,-]','[-,+]','[-,-]'),lty=c(1,1,1,1),col=c('goldenrod3','darkgreen','darkred','black'), text.col=c('goldenrod3','darkgreen','darkred','black'), lwd=c(2,2,2), cex=0.85)
box(lwd=2)
print(fractions)

dev.off()
