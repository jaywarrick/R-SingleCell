#rm(list = ls())

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

source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/logicle.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/plotFACS.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/getSingleStats.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/getDoubleStats.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/drawStats.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R FACS Plotting Functions/drawLogicleAxis.R')
source('/Users/jaywarrick/GoogleDrive/SingleCell/Figures/plotHelperFunctions.R')
library(flowCore)
library(flowViz)
library(MASS)
library(RColorBrewer)
library(foreign)
library(hash)

microscopeFolder <- '/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/Compiled Data Folders/Microscope Data/File - Data Table/'

##### Make M51R DC Movie for ESI #####
setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots')
data <- read.arff(paste(microscopeFolder, 'x1_y1.arff', sep=''))
temp <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
medC <- c(median(temp$G), median(temp$R))
madC <- c(mad(temp$G), mad(temp$R))
threshC <- medC + 3*madC
print(threshC)

for(t in unique(data$Time))
{
    temp <- data.frame(G=subset(data, Time==t & Measurement=='G')$Value, R=subset(data, Time==t & Measurement=='R')$Value)
    pdf(paste('./DC_M51R/Microscope_M51R_',sprintf('%03.f',as.numeric(t)),'.pdf',sep=''),width=1.1*W,height=1.1*H)
    paperParams(1,1)
    dTime <- sprintf('%0.2f',(2 + 5/6) + (1/6)*as.numeric(as.character(t)))
    plotFACS('dummy', temp, xThresh=threshC[1], yThresh=threshC[2], xTransition=threshC[1], yTransition=threshC[2], time=dTime)
    dev.off()
}

setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots/DC_M51R')
cmd1 <- 'convert -verbose -background white -density 300 *.pdf -quality 100 %03d.jpg'
cmd2 <- 'ffmpeg -y -r 12 -i %03d.jpg -vcodec mjpeg -qscale 1 -an output.avi'
system(cmd1)
system(cmd2)

##### Make WT DC Movie for ESI #####
setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots')
data <- read.arff(paste(microscopeFolder, 'x0_y0.arff', sep=''))
temp <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
medC <- c(median(temp$G), median(temp$R))
madC <- c(mad(temp$G), mad(temp$R))
threshC <- medC + 3*madC
print(threshC)

for(t in unique(data$Time))
{
    temp <- data.frame(G=subset(data, Time==t & Measurement=='G')$Value, R=subset(data, Time==t & Measurement=='R')$Value)
    pdf(paste('./DC_WT/Microscope_WT_',sprintf('%03.f',as.numeric(t)),'.pdf',sep=''),width=1.1*W,height=1.1*H)
    paperParams(1,1)
    dTime <- sprintf('%0.2f',(2 + 5/6) + (1/6)*as.numeric(as.character(t)))
    plotFACS('dummy', temp, xThresh=threshC[1], yThresh=threshC[2], xTransition=threshC[1], yTransition=threshC[2], time=dTime)
    dev.off()
}

setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots/DC_WT')
cmd1 <- 'convert -verbose -background white -density 300 *.pdf -quality 100 %03d.jpg'
cmd2 <- 'ffmpeg -y -r 6 -i %03d.jpg -vcodec mjpeg -qscale 1 -an output.avi'
system(cmd1)
system(cmd2)

##### Make WT DC Movie for ESI #####
setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots')
data <- read.arff(paste(microscopeFolder, 'x1_y2.arff', sep=''))
temp <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
medC <- c(median(temp$G), median(temp$R))
madC <- c(mad(temp$G), mad(temp$R))
threshC <- medC + 3*madC
print(threshC)

for(t in unique(data$Time))
{
    temp <- data.frame(G=subset(data, Time==t & Measurement=='G')$Value, R=subset(data, Time==t & Measurement=='R')$Value)
    pdf(paste('./DC_MOCK/Microscope_MOCK_',sprintf('%03.f',as.numeric(t)),'.pdf',sep=''),width=1.1*W,height=1.1*H)
    paperParams(1,1)
    dTime <- sprintf('%0.2f',(2 + 5/6) + (1/6)*as.numeric(as.character(t)))
    plotFACS('dummy', temp, xThresh=threshC[1], yThresh=threshC[2], xTransition=threshC[1], yTransition=threshC[2], time=dTime)
    dev.off()
}

setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots/DC_MOCK')
cmd1 <- 'convert -verbose -background white -density 300 *.pdf -quality 100 %03d.jpg'
cmd2 <- 'ffmpeg -y -r 6 -i %03d.jpg -vcodec mjpeg -qscale 1 -an output.avi'
system(cmd1)
system(cmd2)

##### Make Summary Plots for Manuscript #####
# read in the appropriate microscope data
# We want times 19, 55, and 91 from the microscope data

results <- list()

data <- read.arff(paste(microscopeFolder, 'x0_y0.arff', sep=''))
results$A0 <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
results$A1 <- data.frame(G=subset(data, Time=='19' & Measurement=='G')$Value, R=subset(data, Time=='19' & Measurement=='R')$Value)
results$A2 <- data.frame(G=subset(data, Time=='55' & Measurement=='G')$Value, R=subset(data, Time=='55' & Measurement=='R')$Value)
results$A3 <- data.frame(G=subset(data, Time=='91' & Measurement=='G')$Value, R=subset(data, Time=='91' & Measurement=='R')$Value)

data <- read.arff(paste(microscopeFolder, 'x1_y1.arff', sep=''))
results$C0 <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
results$C1 <- data.frame(G=subset(data, Time=='19' & Measurement=='G')$Value, R=subset(data, Time=='19' & Measurement=='R')$Value)
results$C2 <- data.frame(G=subset(data, Time=='55' & Measurement=='G')$Value, R=subset(data, Time=='55' & Measurement=='R')$Value)
results$C3 <- data.frame(G=subset(data, Time=='91' & Measurement=='G')$Value, R=subset(data, Time=='91' & Measurement=='R')$Value)

data <- read.arff(paste(microscopeFolder, 'x1_y2.arff', sep=''))
results$M0 <- data.frame(G=subset(data, Time=='0' & Measurement=='G')$Value, R=subset(data, Time=='0' & Measurement=='R')$Value)
results$M1 <- data.frame(G=subset(data, Time=='19' & Measurement=='G')$Value, R=subset(data, Time=='19' & Measurement=='R')$Value)
results$M2 <- data.frame(G=subset(data, Time=='55' & Measurement=='G')$Value, R=subset(data, Time=='55' & Measurement=='R')$Value)
results$M3 <- data.frame(G=subset(data, Time=='91' & Measurement=='G')$Value, R=subset(data, Time=='91' & Measurement=='R')$Value)

medA <- c(median(results$A0$G), median(results$A0$R))
madA <- c(mad(results$A0$G), mad(results$A0$R))
threshA <- medA + 3*madA
print(threshA)

medC <- c(median(results$C0$G), median(results$C0$R))
madC <- c(mad(results$C0$G), mad(results$C0$R))
threshC <- medC + 3*madC
print(threshC)

items <- names(results)
for(item in items)
{
    if(length(grep('A',item,fixed=TRUE)) > 0)
    {
        # then it's A data
        
        pdf(paste('WT_',item,'.pdf',sep=''),width=1.1*W,height=1.1*H)
        paperParams(1,1)
        plotFACS('dummy', results[[item]], xThresh=threshA[1], yThresh=threshA[2], xTransition=threshA[1], yTransition=threshA[2])
        dev.off()
    } else if(length(grep('M',item,fixed=TRUE)) > 0)
    {
        # then it's M data
        
        pdf(paste('MOCK_',item,'.pdf',sep=''),width=1.1*W,height=1.1*H)
        paperParams(1,1)
        plotFACS('dummy', results[[item]], xThresh=threshA[1], yThresh=threshA[2], xTransition=threshA[1], yTransition=threshA[2])
        dev.off()
    } else
    {
        # its' C data
        pdf(paste('M51R_',item,'.pdf',sep=''),width=1.1*W,height=1.1*H)
        paperParams(1,1)
        plotFACS('dummy', results[[item]], xThresh=threshC[1], yThresh=threshC[2], xTransition=threshC[1], yTransition=threshC[2])
        dev.off()
    }
}
