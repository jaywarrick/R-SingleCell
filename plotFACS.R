# Use commenting to select the appropriate function header. The first is used for plotting 'FACS' data while the second is used for plotting 'Microscope' data.
# Doing it this way ensures the same visual parameters and labels are used for both types of plots.
# Crossover for FACS = 0.002
# Crossover for Kalin Scope = 0.004
# Crossover for Yin Scope = 0.0007
#plotFACS <- function(title='title', temp, xmin=-100, xmax=200000, ymin=-200, ymax=300000, xThresh=15, yThresh=15, xCross=0, yCross=0.002, xTransition=xThresh, yTransition=yThresh, xLinLogRatio=200, yLinLogRatio=xLinLogRatio, xTicks=c(0,100,1000,10000,100000,1000000), yTicks=xTicks, xTickLabels=expression(0,10^2,10^3,10^4,10^5,10^6), yTickLabels=xTickLabels, time='-1')
#plotFACS <- function(title='title', temp, xmin=-20, xmax=60000, ymin=-30, ymax=80000, xThresh=15, yThresh=15, xCross=0, yCross=0.004, xTransition=xThresh, yTransition=yThresh, xLinLogRatio=50, yLinLogRatio=xLinLogRatio, xTicks=c(0,10,100,1000,10000,100000,1000000), yTicks=xTicks, xTickLabels=expression(0,10^1,10^2,10^3,10^4,10^5,10^6), yTickLabels=xTickLabels, xcol='darkgreen', ycol='darkred', time='-1')
plotFACS <- function(title='title', temp, xmin=-5, xmax=60000, ymin=-30, ymax=100000, xThresh=15, yThresh=15, xCross=0, yCross=0.0007, xTransition=xThresh, yTransition=yThresh, xLinLogRatio=50, yLinLogRatio=xLinLogRatio, xTicks=c(0,10,100,1000,10000,100000,1000000), yTicks=xTicks, xTickLabels=expression(0,10^1,10^2,10^3,10^4,10^5,10^6), yTickLabels=xTickLabels, xcol='darkgreen', ycol='darkred', time='-1')
{   
    #     source('/Users/jaywarrick/Documents/Yin_Lab/Papers/SingleCell/Figures/paperPlot.R')
    #xLinLogRatio=100, yLinLogRatio=xLinLogRatio, xTicks=c(-50,0,15,50,100,1000,10000,100000), yTicks=xTicks)
    
    #     rm(title,temp,xmin,xmax,ymin,ymax,xThresh,yThresh,xCross,yCross,xTransition,yTransition,xLinLogRatio,yLinLogRatio,xTicks,yTicks,xTickLabels,yTickLabels,xcol,ycol)
    #     title<-'title'
    #     xmin <- -20
    #     xmax <- 60000
    #     ymin <- -30
    #     ymax <- 80000
    #     xThresh <- 15
    #     yThresh <- 15
    #     xCross <- 0
    #     yCross <- 0.004
    #     xTransition <- xThresh
    #     yTransition <- yThresh
    #     xLinLogRatio <- 50
    #     yLinLogRatio <- xLinLogRatio
    #     xTicks <- c(0,10,100,1000,10000,100000,1000000)
    #     yTicks <- xTicks
    #     xTickLabels <- expression(0,10^1,10^2,10^3,10^4,10^5,10^6)
    #     yTickLabels <- xTickLabels
    #     xcol <- 'darkgreen'
    #     ycol <- 'darkred'
    
    x <- temp$G.BC.NULL
    y <- temp$R.BC.NULL
    if(is.null(x))
    {
        x <- temp$G
        y <- temp$R
    }
    myStats <- getDoubleStats(x, y, xThresh, xCross, yThresh, yCross)
    print(myStats)
    
    xminAdj <- logicle(xmin, xTransition, xLinLogRatio);
    xmaxAdj <- logicle(xmax, xTransition, xLinLogRatio);
    yminAdj <- logicle(ymin, yTransition, yLinLogRatio);
    ymaxAdj <- logicle(ymax, yTransition, yLinLogRatio);
    print(xminAdj)
    print(xmaxAdj)
    print(yminAdj)
    print(ymaxAdj)
    
    paperParams(1,1)
    par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
    par(mgp = c(2,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
    paperPlot(letter='',xlabel='(Host) GFP Intensity [au]', ylabel='(Virus) RFP Intensity [au]', xlimit=c(xminAdj, xmaxAdj), ylimit=c(yminAdj,ymaxAdj), plotaxes=FALSE, log='xy')
    
    x <- logicle(x, xTransition, xLinLogRatio)
    y <- logicle(y, yTransition, yLinLogRatio)
    
    points(x, y, pch=21, cex=0.35, bg='black')
    
    # Horizontal Thresh (Y)
    tempX <- c(seq(xmin, xTransition, length.out=100), 10^(seq(log10(xTransition), log10(xmax), length.out=100)))
    tempY <- xCross*tempX + yThresh
    tempX <- logicle(tempX, xTransition, xLinLogRatio)
    tempY <- logicle(tempY, yTransition, yLinLogRatio)
    lines(tempX, tempY, col='blue', lwd=1.5)
    
    # Vertical Thresh (X)
    tempY <- c(seq(0.3*ymin, yTransition, length.out=100), 10^(seq(log10(yTransition), log10(0.8*ymax), length.out=100)))
    tempX <- yCross*tempY + xThresh
    tempX <- logicle(tempX, xTransition, xLinLogRatio)
    tempY <- logicle(tempY, yTransition, yLinLogRatio)
    lines(tempX, tempY, col='blue', lwd=1.5)
    
    # Set the axis value
    drawLogicleAxis(xTicks, xTickLabels, 1, xTransition, xLinLogRatio)
    drawLogicleAxis(yTicks, yTickLabels, 2, yTransition, yLinLogRatio)
    
    # Draw stats
    drawStats(myStats, xmin, xmax, ymin, ymax, xTransition, xLinLogRatio, yTransition, yLinLogRatio)
    
    # Draw box
    box(col='black',lwd=1.5)
    
    # Draw time stamp if given
    if(time != '-1')
    {
        mtext(paste('t = ', time, ' [hpi]', sep=''), outer = TRUE, side = 1, adj = 0.05, line = 0, padj = -2.5, cex = 1, col = 'black', font = NA)
    }
    
    return(myStats)
}