drawStats <- function(stats, xmin, xmax, ymin, ymax, xTransition, xLinLogRatio, yTransition, yLinLogRatio)
{

    statsColor = 'blue'
    annotationsize = 1

    n <- 0 # number of decimals to show FOR AVERAGES
    n2 <- 2 # number of decimals to show FOR PERCENTS
    # format(round(x, n), nsmall = n)

    xmin = logicle(xmin, xTransition, xLinLogRatio)# * 5;
    ymin = logicle(ymin, yTransition, yLinLogRatio)# * 1.5;
    xmax = logicle(xmax, xTransition, xLinLogRatio)# * 0.2;
    ymax = logicle(ymax, yTransition, yLinLogRatio)# * 0.66;

    plusplus = NULL
    minusminus = NULL
    plusminus = NULL
    minusplus = NULL

    mean = stats$'XY meanX ++'
    if(is.null(mean) || is.nan(mean))
    {
        plusplus = '';
    }
    else
    {
        plusplus = paste('avg=(', format(round(mean, n), nsmall = n), ',', format(round(stats$'XY meanY ++', n), nsmall = n), ')', sep='')
    }
    mean = stats$'XY meanX --'
    if(is.null(mean) || is.nan(mean))
    {
        minusminus = '';
    }
    else
    {
        minusminus = paste('avg=(', format(round(mean, n), nsmall = n), ',', format(round(stats$'XY meanY --', n), nsmall = n), ')', sep='')
    }
    mean = stats$'XY meanX +-'
    if(is.null(mean) || is.nan(mean))
    {
        plusminus = '';
    }
    else
    {
        plusminus = paste('avg=(', format(round(mean, n), nsmall = n), ',', format(round(stats$'XY meanY +-', n), nsmall = n), ')', sep='')
    }
    mean = stats$'XY meanX -+'
    if(is.null(mean) || is.nan(mean))
    {
        minusplus = '';
    }
    else
    {
        minusplus = paste('avg=(', format(round(mean, n), nsmall = n), ',', format(round(stats$'XY meanY -+', n), nsmall = n), ')', sep='')
    }

    text(xmax, ymax, paste(plusplus, '\n[+,+] ', format(round(stats$"XY % ++", n2), nsmall = n2), '%', sep=''), col = 'goldenrod3', cex = annotationsize, adj = c(1,1))
    text(xmin, ymin, paste(minusminus, '\n[-,-] ', format(round(stats$"XY % --", n2), nsmall = n2), '%', sep=''), col = 'black', cex = annotationsize, adj = c(0,0))
    text(xmax, ymin, paste(plusminus, '\n[+,-] ', format(round(stats$"XY % +-", n2), nsmall = n2), '%', sep=''), col = 'darkgreen', cex = annotationsize, adj = c(1,0))
    text(xmin, ymax, paste(minusplus, '\n[-,+] ', format(round(stats$"XY % -+", n2), nsmall = n2), '%', sep=''), col = 'darkred', cex = annotationsize, adj = c(0,1))

    # Draw the population mean, blue cross-hair
    print(paste('Hello There:', stats$'meanX', ':', stats$'meanY'))
    points(logicle(stats$'meanX',xTransition,xLinLogRatio), logicle(stats$'meanY',yTransition,yLinLogRatio), col='gray70', pch=10, cex=3, lwd=4)

}