paperParams <- function(rows, cols, labsize=1)
{
    par(mfrow = c(rows,cols))
    par(lwd = 1, col='black')
    par(mar = c(3.7,3.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
    par(mgp = c(2.5,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
    par(oma = c(0.1,0.1,0.1,0.1)); # cutoff beyond other measures
    
    par(cex = 1.25)
    par(cex.lab = labsize) # axislabelsize"
    par(cex.axis = 0.75)
}

st <- function(...)
{
    out <- '';
    for(txt in list(...))
    {
        out <- paste(out, as.character(txt), sep='')
    }
    return(out)
}


paperPlot <- function(xlimit=c(0.01,1), ylimit=c(0.01,1), xlabel='x', ylabel='y', xcol='darkgreen', ycol='darkred', plotaxes=TRUE, log='', xticks=NULL, yticks=NULL, letter='')
{
    
    
    #     xlimit=c(min(lx),max(lx))
    #     ylimit=c(min(ly),max(ly))
    #     xlabel=expression(paste('[Virus] ', alpha, ' (1/h)'))
    #     ylabel='[Virus] Max Intensity'
    #     xcol=xcol
    #     ycol=ycol
    #     log='y'
    #     plot(x=c(),y=c(), xlim=xlimit, ylim=ylimit, xlab=xlabel, ylab=ylabel, axes=TRUE, log=log)
    #     
    plot(c(),c(), xlim=xlimit, ylim=ylimit, xlab='', ylab='', axes=FALSE, log=log)
    title(ylab=ylabel, col.lab=ycol)
    title(xlab=xlabel, col.lab=xcol)
    if(plotaxes)
    {
        
        if(is.null(yticks))
        {
            axis(2, las=1)
        }
        else if(is.character(yticks))
        {
            axis(2, las=1, at=as.numeric(yticks), labels=yticks)
        }
        else if(is.expression(yticks))
        {
            axis(2, las=1, at=unlist(lapply(eval(as.list(yticks)), 'eval')), labels=yticks)
        }
        
        if(is.null(xticks))
        {
            axis(1, las=1)
        }
        else if(is.character(xticks))
        {
            axis(1, las=1, at=as.numeric(xticks), labels=xticks)
        }
        else if(is.expression(xticks))
        {
            axis(1, las=1, at=unlist(lapply(eval(as.list(xticks)), 'eval')), labels=xticks)
        }
    }
    
    box(col='black',lwd=2)
    if(letter != '')
    {
        mtext(side=3,adj=0,text=paste('(',letter,')',sep=''),padj=1.8, outer=TRUE, cex=1.5)  
    }
}

paperHist2 <- function(letter='a)', histogram.R, histogram.G, xlabel='bins', ylabel='Density', labsize=1.5, xlim, ylim)
{
    par(mar = c(3.7,4,0.8,0.8)); # beyond plot frame [b,l,t,r]
    plot(getHistStuff(histogram.R), xlab=xlabel, ylab=ylabel, xlim=xlim, ylim=c(min(ylim), 1.2*max(ylim)), xaxs='i', yaxs='i', main=NULL, axes=FALSE, type='l', col=NA)
    box(col='black',lwd=2, bty='l')
    polygon(getHistStuff(histogram.R), col=rgb(1,0,0,0.4))
    polygon(getHistStuff(histogram.G), col=rgb(0,1,0,0.2))
    axis(1,las=1)
    axis(2,las=1)
    box(col='black',lwd=2, bty='l')
    mtext(side=3,adj=0.04,text=paste('(',letter,')',sep=''),padj=1.8, outer=TRUE, cex=1.3)
}

getHistStuff <- function(histogram)
{
    dx <- histogram$mids[2] - histogram$mids[1]
    xL <- histogram$mids - 0.5*dx
    xR <- histogram$mids + 0.5*dx
    x <- xL[1];
    y <- 0;
    for(i in 1:length(xL))
    {
        x <- c(x, xL[i], xR[i]) 
        y <- c(y, histogram$density[i], histogram$density[i]) 
    }
    x <- c(x, xR[length(xR)])
    y <- c(y, 0)
    return(data.frame(x=x, y=y))
}

paperHist <- function(letter='a)',histogram, xlabel='bins', ylabel='Density', labsize=1.5, xlim, ylim, xcol='red', bar.col='grey')
{
    #     plot(c(), c(), frame.plot=FALSE, xlab=xlabel, ylab='Density', cex.lab=labsize, main=NULL, yaxs='i', ylim=c(0,max(histogram$density)*1.04), xlim=c(min(histogram$breaks), max(histogram$breaks)))
    #     lines(histogram) 
    par(mar = c(3.7,4,0.8,0.8)); # beyond plot frame [b,l,t,r]
    myHist(histogram, main=NULL, axes=FALSE, xlab=xlabel, ylab=ylabel, xlim=xlim, ylim=ylim, freq=FALSE, col=bar.col, col.lab=xcol)
    axis(1,las=1)
    axis(2,las=1)
    box(col='black',lwd=2, bty='l')
    mtext(side=3,adj=0,text=paste('(',letter,')',sep=''),padj=1.8, outer=TRUE, cex=1.3)
}

myHist <- function (histogram, freq = equidist, density = NULL, angle = 45, col = NULL, 
                    border = par("fg"), lty = NULL, main = paste("Histogram of", 
                                                                 histogram$xname), xlim = range(histogram$breaks), ylim = NULL, xlab = histogram$xname, 
                    ylab, axes = TRUE, labels = FALSE, add = FALSE, width=1.0, offset=(1.0-width)/2, ...) 
{
    y <- histogram$counts
    if(!freq)
    {
        y <- histogram$density
    }
    nB <- length(histogram$breaks)
    if (is.null(ylim)) 
        ylim <- range(y, 0)
    if (missing(ylab)) 
        ylab <- if (!freq) 
            "Density"
    else "Frequency"
    plot.new()
    plot.window(xlim, ylim, "", xaxs='i', yaxs='i')
    title(main = main, xlab = xlab, ylab = ylab, ...)
    if (axes) {
        axis(1, ...)
        axis(2, ...)
    }
    
    if (width != 1.0 || offset != 0) {
        # Calculates the width of each bar in the histogram
        delta.breaks <- histogram$breaks[-1] - histogram$breaks[-nB];
        x.offset <- offset * delta.breaks;
        x.width <- width * delta.breaks;
        x <- histogram$breaks[-nB]+x.offset;
        rect(x, 0, x+x.width, y, col=col, border=border, angle = angle, density = density, lty=lty);
    } else {
        rect(histogram$breaks[-nB], 0, histogram$breaks[-1], y, col = col, border = border, 
             angle = angle, density = density, lty = lty)
    }
    
    if ((logl <- is.logical(labels) && labels) || is.character(labels)) 
        text(histogram$mids, y, labels = if (logl) {
            if (freq) 
                histogram$counts
            else round(histogram$density, 3)
        }
        else labels, adj = c(0.5, -0.5))
    invisible()
}

DefinePWFct <- function(x1, x2)
{
    # Define a piecewise linear model
    pNames.pw <- c('L','R');
    func <- paste('function(x, par)
    {
                      y1 <- par[,1]
                      y2 <- par[,2]
                      y <- x
                      
                      left <- which(x <=', x1, ');
                      y[left] <- y1[left];
                      
                      right <- which(x >=', x2, ');
                      y[right] <- y2[right];
                      
                      m <- (y2-y1)/(', x2, '-', x1, ')
                      b <- y1 - m*x1
                      
                      middle <- which(x >', x1, ' & x <',  x2, ')
                      y[middle] <- m[middle]*x[middle] + b[middle]
                      
                      y
    }',sep = '')
    
    fct.pw <- eval(parse(text=func))
    
    func <- paste('function(x, par)
    {
                      y1 <- par[1]
                      y2 <- par[2]
                      y <- x
                      
                      left <- which(x <=', x1, ');
                      y[left] <- y1;
                      
                      right <- which(x >=', x2, ');
                      y[right] <- y2;
                      
                      m <- (y2-y1)/(', x2, '-', x1, ')
                      b <- y1 - m*x1
                      
                      middle <- which(x >', x1, ' & x <',  x2, ')
                      y[middle] <- m*x[middle] + b
                      
                      y
    }',sep = '')
    
    fct.pw2 <- eval(parse(text=func))
    
    ssfct.pw <- function(data) 
    { 
        return(c(-1,1)) 
    }
    fctDef.pw <- list(name='pwlin',fct=fct.pw, fct2=fct.pw2, ssfct=ssfct.pw, names=pNames.pw);
    return(list(a=pNames.pw, b=fct.pw, b2=fct.pw2, c=ssfct.pw, d=fctDef.pw))
}

daFunc <- function(x, par)
{
    myy1 <- par[,1]
    myy2 <- par[,2]
    myx1 <- par[,3]
    myx2 <- par[,4]
    
    #     x <- pop$Time.Delta.Delay
    myy <- x
    
    left <- which(x <= myx1);
    myy[left] <- myy1[left];
    
    right <- which(x >= myx2);
    myy[right] <- myy2[right];
    
    m <- (myy2-myy1)/(myx2 - myx1)
    b <- myy1 - m*myx1
    
    middle <- which(x > myx1 & x < myx2)
    myy[middle] <- m[middle]*x[middle] + b[middle]
    
    myy
}

daFunc2 <- function(x, par)
{
    myy1 <- par[1]
    myy2 <- par[2]
    myx1 <- par[3]
    myx2 <- par[4]
    
    #     x <- pop$Time.Delta.Delay
    myy <- x
    
    left <- which(x <= myx1);
    myy[left] <- myy1;
    
    right <- which(x >= myx2);
    myy[right] <- myy2;
    
    m <- (myy2-myy1)/(myx2 - myx1)
    b <- myy1 - m*myx1
    
    middle <- which(x > myx1 & x < myx2)
    myy[middle] <- m*x[middle] + b
    
    myy
}

startDaFunc <- function(data) 
{ 
    return(c(-1, 1, -8, 6)) 
    #     return(c(-4.711, -0.9105, 6, 1.181)) 
}

daFuncDef <- list(name='pwlin',fct=daFunc, ssfct=startDaFunc, names=c('Y_1','Y_2','X_1','X_2'));

getDensityColors <- function(x, y)
{
    #     library(hexbin)
    xunique <- unique(x)
    yunique <- unique(y)
    xync <- data.frame(x=c(), y=c(), n=c(), c=c())
    for(xu in xunique)
    {
        for(yu in yunique)
        {
            i <- sum(x==xu & y==yu)
            xync <- rbind(xync, data.frame(x=xu, y=yu, n=i))
        }
    }
    xync <- subset(xync, n > 0)
    theColors <- rev(rainbow(max(xync$n), end=0.6))
    #     theColors <- plinrain(5, beg=50, end=200)
    #     plot(c(), c(), xlim=c(0, length(theColors) + 1), ylim=c(0, length(theColors) + 1))
    #     for(i in 1:length(theColors))
    #     {
    #         points(i, i, pch=21, col=rgb(0,0,0,0), bg=theColors[i])
    #     }
    xync$c <- theColors[xync$n]
    return(list(colors=theColors, data=xync))
}