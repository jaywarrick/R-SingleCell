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

FACSFolder <- '/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/Compiled Data Folders/FCS Data/Compiled/'
setwd('/Users/jaywarrick/GoogleDrive/SingleCell/FACS Data/R Plots/')

thresh <- list()
med <- list()
files <- c('WT_6_B.fcs', 'WT_12_C.fcs', 'WT_S_A.fcs', 'M51R_6_A.fcs', 'M51R_12_C.fcs', 'M51R_S_B.fcs', 'Mock_6_C.fcs', 'Mock_12_C.fcs', 'Mock_S_B.fcs')
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
    pdf(paste(f,'.pdf',sep=''),width=1.1*W,height=1.1*H)
    plotFACS(f, temp, xThresh=gateX, yThresh=gateY, xTransition=gateX, yTransition=gateY)
    dev.off()
}


