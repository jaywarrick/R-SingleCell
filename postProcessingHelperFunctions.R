library(drc)
library(zoo)
library(plyr)

getNoise <- function(temp, X='R', n=3)
{
    return(min(temp[temp[,X]>0,X])/3)
    #     times <- unique(temp$time)
    #     times <- sort(times)[1:n]
    #     return(sd(subset(temp, time %in% times)[,X]))
}

getDeath <- function(temp, X='R', window=5, thresh=0.6)
{
    require(zoo)
    roll <- rollmean(temp[,X], k=window, na.pad=TRUE, align='right')
    deaths <- which(temp[,X][2:length(temp$time)] < thresh*roll[1:(length(temp$time)-1)])
    if(length(deaths) > 0)
    {
        # if death is detected 
        return(list(t=temp$time[deaths[1]], int=temp[,X][deaths[1]]))
    }
    else
    {
        return(NA)
    }
}

getTimeDelayIndex <- function(temp, X='R')
{
    temp2 <- which(temp[,X] == 0)
    
    if(is.na(temp2[1]))
    {
        # if no values are 0 the first index is the TimeDelayIndex
        return(1)
    }
    
    # find the end of the first set of consecutive 0's and return the next index (i.e., i)
    for(i in 1:(length(temp2)))
    {
        if(i != temp2[i])
        {
            return(i)
        }
    }
    
    # if we reach here then all the 0's were consecutive
    if(length(temp2) == length(temp[,X]))
    {
        # Then the whole trajectory is 0 and there is no TimeDelayIndex
        return(NA)
    }
    else
    {
        # return the index after all the consecutive 0's
        return(length(temp2) + 1)
    }
}

getPeak <- function(temp, X='R', thresh=0.85)
{
    maxX <- max(temp[,X])
    index <- which(temp[,X] >= thresh*maxX)[1]
    return(list(t=temp$time[index], int=temp[,X][index]))
}

getWeightedR2 <- function(y, model)
{
    r <- residuals(model)
    f <- fitted(model)
    w <- weights(model)
    return(getWeightedR2_Vectors(y, r, f, w))
}

getWeightedR2_Vectors <- function(y, r, f, w)
{
    SSr <- sum(w*r^2);
    SSt <- sum(w*(y-mean(y))^2)
    return(1-(SSr/SSt))
}

getLMFit <- function(tempR, tempG)
{    
    LMFit.R <- lm(log(tempR$R)~tempR$time)
    LMFit.R2.R <- summary(LMFit.R)$r.squared
    LMFit.Alpha.R <- coef(LMFit.R)[2]
    LMFit.Tau.R <- -1*coef(LMFit.R)[1]/LMFit.Alpha.R
    LMFit.Flag.Quality.R <- (LMFit.R2.R >= 0.7 & LMFit.Alpha.R > 0)
    
    LMFit.G <- NA
    LMFit.R2.G <- NA
    LMFit.Alpha.G <- NA
    LMFit.Tau.G <- NA
    LMFit.Flag.Quality.G <- NA
    if(!is.null(tempG))
    {
        LMFit.G <- lm(log(tempG$G)~tempG$time)
        LMFit.R2.G <- summary(LMFit.G)$r.squared
        LMFit.Alpha.G <- coef(LMFit.G)[2]
        LMFit.Tau.G <- -1*coef(LMFit.G)[1]/LMFit.Alpha.G
        LMFit.Flag.Quality.G <- (LMFit.R2.G >= 0.7 & LMFit.Alpha.G > 0)
    }
    
    return(list(LMFit.R=LMFit.R, LMFit.G=LMFit.G, LMFit.Alpha.R=LMFit.Alpha.R, LMFit.Alpha.G=LMFit.Alpha.G, LMFit.Tau.R=LMFit.Tau.R, LMFit.Tau.G=LMFit.Tau.G, LMFit.R2.R=LMFit.R2.R, LMFit.R2.G=LMFit.R2.G, LMFit.Flag.Quality.R=LMFit.Flag.Quality.R, LMFit.Flag.Quality.G=LMFit.Flag.Quality.G))
}

plotID <- function(p, temp, R, G)
{
    #     print('plotting');
    print(p$ID)

    tmin <- 0
    tmax <- max(temp$time)
    ymin <- 1
    ymax <- max(max(temp$R),max(temp$G))
    
    plot(temp$time, temp$R, col='red', pch=21, bg='red', type='l', log='y', xlim=c(tmin,tmax), ylim=c(ymin,ymax), xlab='Time', ylab='Intensity', main=id, cex.lab=0.7, cex.axis=0.7)
    points(R$time, R$R, col='red', pch=21, bg='red')
    lines(R$time, R$fitted, col='black')
    points(p$Time.Peak.R, p$Int.Peak.R, col='black', pch=21, bg='black')

    if(!is.na(G)) # then we should have tried to fit green as well
    {
        lines(temp$time, temp$G, col='green', pch=21, bg='green')
        points(G$time, G$G, col='green', pch=21, bg='green')
        lines(G$time, G$fitted, col='black')
        points(p$Time.Peak.G, p$Int.Peak.G, col='black', pch=21, bg='black')
    }
    
    if(!is.na(p$Time.Death))
    {
        points(p$Time.Death, p$Int.Death)
    }

    box(lwd=2)
    
    f <- paste('/Users/jaywarrick/GoogleDrive/SingleCell/DRCFitting/',as.character(id),'_R.pdf',sep='')
    dev.copy2pdf(file=f, width=5, height=5)
    
    return(params)
}

getPar <- function(id, par, model)
{
    label <- paste(par, ':', id, sep='')
    return(coef(model)[label])
}

# Define various equations to fit the trajectory data with
# Define y = e^(alpha*(x-tau))
fctA <- function(x, par, BG)
{
    exp(par[, 1] *(x - par[, 2]))
}
ssfctA <- function(data)
{
    # names=c('alpha','tau')
    #     maxy <- max(data$y)[1]
    #     t <- data$x[data$y==maxy]
    #     print(data)
    # Given, the exponential nature of the response, return a simple exponential through the max point because that will have the greatest impact on the least squares residual
    #     return(c(log(maxy)/t,3.5))
    modelo <- lm(log(data$y)~data$x)
    #     print(summary(modelo))
    alpha <- coef(modelo)[2]
    tau <- -1*coef(modelo)[1]/alpha
    return(c(alpha, tau))
}

getDualParams <- function(id, temp, Int.Max.R, Int.Max.G) # for RFP+ GFP+ Microwells
{
    #     id <- IDs[10]
    #     Int.Max.R <- max(temp$R)
    #     Int.Max.G <- max(temp$G)
    Time.Max.R <- temp$time[temp$R==Int.Max.R][1]
    Time.Max.G <- temp$time[temp$G==Int.Max.G][1]
    print(id)
    
    # Find some critical times and values for analysis
    Time.Delay.Index.R <- getTimeDelayIndex(temp, 'R')
    Time.Delay.Index.G <- getTimeDelayIndex(temp, 'G')
    Time.Delay.R <- temp$time[Time.Delay.Index.R]
    Time.Delay.G <- temp$time[Time.Delay.Index.G]
    duh <- getPeak(temp, 'R')
    Time.Peak.R <- duh$t
    Int.Peak.R <- duh$int
    duh <- getPeak(temp, 'G')
    Time.Peak.G <- duh$t
    Int.Peak.G <- duh$int
    duh <- getDeath(temp, 'R')
    if(is.null(duh)[1] | is.na(duh)[1])
    {
        Time.Death <- NA
        Int.Death <- NA
    }else
    {
        Time.Death <- duh$t
        Int.Death <- duh$int            
    }
    
    # Save a version of the data with only the first four "detectable" datapoints
    tempR <- subset(temp, time %in% temp$time[Time.Delay.Index.R:(Time.Delay.Index.R + 3)])
    tempG <- subset(temp, time %in% temp$time[Time.Delay.Index.G:(Time.Delay.Index.G + 3)])
    
    # Log the data and fit it with a line using the lm function to get an estimate of alpha and tau
    LMFit <- list()
    LMFit <- getLMFit(tempR, tempG)
    
    # Save the fit information into a table called params
    newParams <- data.frame(ID=id, x=tempR$x[1], y=tempR$y[1], Device=tempR$Device[1], ImRow=tempR$ImRow[1], ImCol=tempR$ImCol[1], ROI=tempR$ROI[1], Time.Max.R=Time.Max.R, Time.Max.G=Time.Max.G, Time.Delay.R=Time.Delay.R, Time.Delay.G=Time.Delay.G, Time.Peak.R=Time.Peak.R, Time.Peak.G=Time.Peak.G, Time.Death=Time.Death, Int.Max.R=Int.Max.R, Int.Max.G=Int.Max.G, Int.Peak.R=Int.Peak.R, Int.Peak.G=Int.Peak.G, Int.Death=Int.Death, LMFit[3:length(LMFit)])
    return(list(newDataR=tempR, newDataG=tempG, newParams=newParams))
}

getSingleParams <- function(id, temp, Int.Max.R) # for RFP+ GFP- Microwells
{
    Time.Max.R <- temp$time[temp$R==Int.Max.R][1]
    Time.Max.G <- 0
    print(id)
    
    # Find some critical times and values for analysis
    Time.Delay.Index.R <- getTimeDelayIndex(temp, 'R')
    Time.Delay.Index.G <- NA
    Time.Delay.R <- temp$time[Time.Delay.Index.R]
    Time.Delay.G <- NA
    
    duh <- getPeak(temp, 'R')
    Time.Peak.R <- duh$t
    Int.Peak.R <- duh$int
    
    Time.Peak.G <- NA
    Int.Peak.G <- NA
    
    duh <- getDeath(temp, 'R')
    if(is.null(duh)[1] | is.na(duh)[1])
    {
        Time.Death <- NA
        Int.Death <- NA
    }else
    {
        Time.Death <- duh$t
        Int.Death <- duh$int            
    }
    
    # Save a version of the data with only the first four "detectable" datapoints
    tempR <- subset(temp, time %in% temp$time[Time.Delay.Index.R:(Time.Delay.Index.R + 3)])
    
    # Log the data and fit it with a line using the lm function to get an estimate of alpha and tau
    LMFit <- list()
    LMFit <- getLMFit(tempR, NULL)
    
    # Save the fit information into a table called params
    #     print(print('inner'))
    newParams = data.frame(ID=id, x=tempR$x[1], y=tempR$y[1], Device=tempR$Device[1], ImRow=tempR$ImRow[1], ImCol=tempR$ImCol[1], ROI=tempR$ROI[1], Time.Max.R=Time.Max.R, Time.Max.G=Time.Max.G, Time.Delay.R=Time.Delay.R, Time.Delay.G=Time.Delay.G, Time.Peak.R=Time.Peak.R, Time.Peak.G=Time.Peak.G, Time.Death=Time.Death, Int.Max.R=Int.Max.R, Int.Max.G=Int.Max.G, Int.Peak.R=Int.Peak.R, Int.Peak.G=Int.Peak.G, Int.Death=Int.Death, LMFit[3:length(LMFit)])
    #     print(newParams)
    return(list(newDataR=tempR, newDataG=NA, newParams=newParams))
}

siegel.tukey <- function(x, y, id.col = FALSE, adjust.median = F, 
                         rnd = -1, alternative = "two.sided", mu = 0, paired = FALSE, 
                         exact = FALSE, correct = TRUE, conf.int = FALSE, conf.level = 0.95) {
    ###### published on:
    #   http://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/
    ## Main author of the function:  Daniel Malter
    
    # x: a vector of data
    
    # y: Group indicator (if id.col=TRUE); data of the second
    #   group (if
    # id.col=FALSE). If y is the group indicator it MUST take 0
    #   or 1 to indicate
    # the groups, and x must contain the data for both groups.
    
    # id.col: If TRUE (default), then x is the data column and y
    #   is the ID column,
    # indicating the groups. If FALSE, x and y are both data
    #   columns. id.col must
    # be FALSE only if both data columns are of the same length.
    
    # adjust.median: Should between-group differences in medians
    #   be leveled before
    # performing the test? In certain cases, the Siegel-Tukey
    #   test is susceptible
    # to median differences and may indicate significant
    #   differences in
    # variability that, in reality, stem from differences in
    #   medians.
    
    # rnd: Should the data be rounded and, if so, to which
    #   decimal? The default
    # (-1) uses the data as is. Otherwise, rnd must be a
    #   non-negative integer.
    # Typically, this option is not needed. However,
    #   occasionally, differences in
    # the precision with which certain functions return values
    #   cause the merging
    # of two data frames to fail within the siegel.tukey
    #   function. Only then
    # rounding is necessary. This operation should not be
    #   performed if it affects
    # the ranks of observations.
    
    #  arguments passed on to the Wilcoxon test. See
    #   ?wilcox.test
    
    # Value: Among other output, the function returns the data,
    #   the Siegel-Tukey
    # ranks, the associated Wilcoxons W and the p-value for a
    #   Wilcoxon test on
    # tie-adjusted Siegel-Tukey ranks (i.e., it performs and
    #   returns a
    # Siegel-Tukey test). If significant, the group with the
    #   smaller rank sum has
    # greater variability.
    
    # References: Sidney Siegel and John Wilder Tukey (1960) A
    #   nonparametric sum
    # of ranks procedure for relative spread in unpaired
    #   samples. Journal of the
    # American Statistical Association. See also, David J.
    #   Sheskin (2004)
    # Handbook of parametric and nonparametric statistical
    #   procedures. 3rd
    # edition. Chapman and Hall/CRC. Boca Raton, FL.
    
    # Notes: The Siegel-Tukey test has relatively low power and
    #   may, under certain
    # conditions, indicate significance due to differences in
    #   medians rather than
    # differences in variabilities (consider using the argument
    #   adjust.median).
    
    # Output (in this order)
    
    # 1. Group medians (after median adjustment if specified)
    # 2. Wilcoxon-test for between-group differences in medians
    #   (after the median
    # adjustment if specified)
    # 3. Data, group membership, and the Siegel-Tukey ranks
    # 4. Mean Siegel-Tukey rank by group (smaller values indicate
    #   greater
    # variability)
    # 5. Siegel-Tukey test (Wilcoxon test on tie-adjusted
    #   Siegel-Tukey ranks)
    
    is.formula <- function(x) class(x) == "formula"
    
    if (is.formula(x)) {
        y <- do.call(c, list(as.name(all.vars(x)[2])), envir = parent.frame(2))
        x <- do.call(c, list(as.name(all.vars(x)[1])), envir = parent.frame(2))  # I am using parent.frame(2) since if the name of the variable in the equation is 'x', then we will mistakenly get the function in here instead of the vector.
        id.col <- TRUE
        # print(x)
        # print(ls.str())
        # data=data.frame(c(x,y),rep(c(0,1),c(length(x),length(y))))
        # print(data)
    }
    
    if (id.col == FALSE) {
        data = data.frame(c(x, y), rep(c(0, 1), c(length(x), length(y))))
    } else {
        data = data.frame(x, y)
    }
    names(data) = c("x", "y")
    data = data[order(data$x), ]
    if (rnd > -1) {
        data$x = round(data$x, rnd)
    }
    
    if (adjust.median == T) {
        cat("\n", "Adjusting medians...", "\n", sep = "")
        data$x[data$y == 0] = data$x[data$y == 0] - (median(data$x[data$y == 
                                                                       0]))
        data$x[data$y == 1] = data$x[data$y == 1] - (median(data$x[data$y == 
                                                                       1]))
    }
    cat("\n", "Median of group 1 = ", median(data$x[data$y == 0]), 
        "\n", sep = "")
    cat("Median of group 2 = ", median(data$x[data$y == 1]), "\n", 
        "\n", sep = "")
    cat("Testing median differences...", "\n")
    print(wilcox.test(data$x[data$y == 0], data$x[data$y == 1]))
    
    # The following must be done for the case when id.col==F
    x <- data$x
    y <- data$y
    
    cat("Performing Siegel-Tukey rank transformation...", "\n", 
        "\n")
    
    
    
    sort.x <- sort(data$x)
    sort.id <- data$y[order(data$x)]
    
    data.matrix <- data.frame(sort.x, sort.id)
    
    base1 <- c(1, 4)
    iterator1 <- matrix(seq(from = 1, to = length(x), by = 4)) - 
        1
    rank1 <- apply(iterator1, 1, function(x) x + base1)
    
    iterator2 <- matrix(seq(from = 2, to = length(x), by = 4))
    base2 <- c(0, 1)
    rank2 <- apply(iterator2, 1, function(x) x + base2)
    
    #print(rank1)
    #print(rank2)
    
    if (length(rank1) == length(rank2)) {
        rank <- c(rank1[1:floor(length(x)/2)], rev(rank2[1:ceiling(length(x)/2)]))
    } else {
        rank <- c(rank1[1:ceiling(length(x)/2)], rev(rank2[1:floor(length(x)/2)]))
    }
    
    
    unique.ranks <- tapply(rank, sort.x, mean)
    unique.x <- as.numeric(as.character(names(unique.ranks)))
    
    rank.matrix <- data.frame(unique.x, unique.ranks)
    
    ST.matrix <- merge(data.matrix, rank.matrix, by.x = "sort.x", 
                       by.y = "unique.x")
    
    print(ST.matrix)
    
    cat("\n", "Performing Siegel-Tukey test...", "\n", sep = "")
    
    ranks0 <- ST.matrix$unique.ranks[ST.matrix$sort.id == 0]
    ranks1 <- ST.matrix$unique.ranks[ST.matrix$sort.id == 1]
    
    cat("\n", "Mean rank of group 0: ", mean(ranks0), "\n", sep = "")
    cat("Mean rank of group 1: ", mean(ranks1), "\n", sep = "")
    
    print(wilcox.test(ranks0, ranks1, alternative = alternative, 
                      mu = mu, paired = paired, exact = exact, correct = correct, 
                      conf.int = conf.int, conf.level = conf.level))
    
    return(list(ranks0=ranks0, ranks1=ranks1))
}

y <- rnorm(100)
hist(y)
x <- seq_along(y)
y2 <- rep(y, each=2)
y2 <- y2[-length(y2)]
x2 <- rep(x, each=2)[-1]
x3 <- c(min(x2), x2, max(x2))
y3 <- c(0, y2, 0)

# because polygon() is dumb and wants a pre-existing plot
plot(x, y, ylim=c(0, max(y)), type="n")

polygon(x3, y3, border=NA, col="grey")
lines(x2, y2)


library(ggplot2)
library(Hmisc) # cut2
y1 <- rnorm(100)
y2 <- rnorm(100, mean=1)
breakers <- cut(c(y1,y2), 15)
duh <- hist(c(y1,y2), cuts=breakers)
barplot(table(breaks), space=0, border=NA, col=rgb(0,0,0,1))
duh <- data.frame(y1=y1, y2=y2)
ggplot(duh, aes(x=duh$y1)) + 
    geom_histogram(data=data.frame(y1=duh$y1), fill = "red", alpha = 0.2) + 
    geom_histogram(data=data.frame(y2=duh$y2), fill = "blue", alpha = 0.2)

qplot(y1,data=duh,geom="histogram",fill='red') + qplot(y1,data=duh,geom="histogram",fill='green')
ggplot(guh, aes(length, fill = veg)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')


library(lattice)
types.plain <- c("p", "l", "o", "r", "g", "s", "S", "h", "a", "smooth")
types.horiz <- c("s", "S", "h", "a", "smooth")
horiz <- rep(c(FALSE, TRUE), c(length(types.plain), length(types.horiz)))

types <- c(types.plain, types.horiz)

x <- sample(seq(-10, 10, length.out = 15), 30, TRUE)
y <- x + 0.25 * (x + 1)^2 + rnorm(length(x), sd = 5)

xyplot(y ~ x | gl(1, length(types)),
       xlab = "type", 
       ylab = list(c("horizontal=TRUE", "horizontal=FALSE"), y = c(1/6, 4/6)),
       as.table = TRUE, layout = c(5, 3),
       between = list(y = c(0, 1)),
       strip = function(...) {
           panel.fill(trellis.par.get("strip.background")$col[1])
           type <- types[panel.number()]
           grid::grid.text(label = sprintf('"%s"', type), 
                           x = 0.5, y = 0.5)
           grid::grid.rect()
       },
       scales = list(alternating = c(0, 2), tck = c(0, 0.7), draw = FALSE),
       par.settings = 
           list(layout.widths = list(strip.left = c(1, 0, 0, 0, 0))),
       panel = function(...) {
           type <- types[panel.number()]
           horizontal <- horiz[panel.number()]
           panel.xyplot(..., 
                        type = type,
                        horizontal = horizontal)
       })[rep(1, length(types))]





# getTestData <- function(dataPar, n)
# {
#     library(data.table)
#     data <- data.frame(x=1:7)
#     data$y <- exp(dataPar[1]*(data$x - dataPar[2]))
#     data$weights <- 1/((0.059^2)*(data$y + 350))
#     data <- data.frame(data, ID=gl(n,nrow(data)))
#     data <- data.table(data)
#     data$random <- 1
#     #     browser()
#     data[,random:=rnorm(1,mean=1, sd=0.1),by=ID]
#     data$y <- data$y*data$random + 0.1*data$x
#     data <- data.frame(data)
#     return(data)
# }
# 
# testNLFit <- function(data, tol)
# {
#     
#     
#     mod <- drm(y~x,data=data,curveid=ID,pmodels=data.frame(ID,ID),fct=globalFctDefA, weights=data$weights, na.action=na.omit, separate=FALSE, lowerl=c(0,0), control=drmc(relTol=tol))
#     #     print(summary(mod))
#     data$fitted <- fitted(mod)
#     data$residuals <- residuals(mod)
#     data <- data.table(data)
#     plot(c(), c(), type='p', log='y', ylim=range(data$y), xlim=range(data$x))
#     duh <- function(x,y,fitted)
#     {
#         points(x, y, pch=21, cex=0.5)
#         lines(x, fitted, col='red')
#     }
#     data[,duh(x,y,fitted),by=ID]
#     
#     sse <- data[,sum(weights*residuals^2), by=ID]
#     return(list(data=data.frame(data), sse=sse, mod=mod))
# }
# data <- getTestData(c(1,5), 50)
# results <- testNLFit(data, tol=1e-5)
# # results <- testNLFit(c(1,5), 50)
# mean(results$sse$V1)