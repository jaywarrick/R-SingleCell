getDoubleStats <- function(x=1:10, y=1:10, xThresh=5, xCross=0, yThresh=5, yCross=0)
{
    ret <- getSingleStats(TRUE, x, y, xThresh, xCross);
    ret$values <- getSingleStats(FALSE, x, y, yThresh, yCross);

    total <- ret$'n'

    # Define subpopulations
    plusplus   <- ((x > (y*yCross + xThresh)) & (y > (x*xCross + yThresh)))
    minusminus <- ((x <= (y*yCross + xThresh)) & (y <= (x*xCross + yThresh)))
    plusminus  <- ((x > (y*yCross + xThresh)) & (y <= (x*xCross + yThresh)))
    minusplus  <- ((x <= (y*yCross + xThresh)) & (y > (x*xCross + yThresh)))

    # Calculate stats for each subpopulation
    # ++
    ret$'XY % ++' <- 100 * (sum(plusplus) / (total))
    ret$'XY n ++' <- sum(plusplus)
    ret$'XY meanX ++' <- mean(x[plusplus])
    ret$'XY meanY ++' <- mean(y[plusplus])
    ret$'XY sdX ++' <- sd(x[plusplus])
    ret$'XY sdY ++' <- sd(y[plusplus])

    # --
    ret$'XY % --' <- 100 * (sum(minusminus) / (total))
    ret$'XY n --' <- sum(minusminus)
    ret$'XY meanX --' <- mean(x[minusminus])
    ret$'XY meanY --' <- mean(y[minusminus])
    ret$'XY sdX --' <- sd(x[minusminus])
    ret$'XY sdY --' <- sd(y[minusminus])

    # -+
    ret$'XY % -+' <- 100 * (sum(minusplus) / (total))
    ret$'XY n -+' <- sum(minusplus)
    ret$'XY meanX -+' <- mean(x[minusplus])
    ret$'XY meanY -+' <- mean(y[minusplus])
    ret$'XY sdX -+' <- sd(x[minusplus])
    ret$'XY sdY -+' <- sd(y[minusplus])

    # +-
    ret$'XY % +-' <- 100 * (sum(plusminus) / (total))
    ret$'XY n +-' <- sum(plusminus)
    ret$'XY meanX +-' <- mean(x[plusminus])
    ret$'XY meanY +-' <- mean(y[plusminus])
    ret$'XY sdX +-' <- sd(x[plusminus])
    ret$'XY sdY +-' <- sd(y[plusminus])

    return(ret);
}