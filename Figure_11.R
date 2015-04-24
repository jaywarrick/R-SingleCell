# Figure 2-6: (Include Red lysed dots)

my.prop.test <- function(xThresh=0, x1, x2, y1, y2)
{
    lo.1 <- which(x1 < xThresh)
    hi.1 <- which(x1 >= xThresh)
    lo.2 <- which(x2 < xThresh)
    hi.2 <- which(x2 >= xThresh)
    
    n1.lo <- length(lo.1)
    n1.hi <- length(hi.1)
    n2.lo <- length(lo.2)
    n2.hi <- length(hi.2)

    mat <- as.matrix(data.frame(intact=c(n1.lo,n1.hi), lysed=c(n2.lo,n2.hi)))
    print(mat)
    ret <- prop.test(mat)
    return(list(matrix=mat, test=ret))
}

######### Plot A ##########

# a) x=Time.Delta.Delay y=Int.Max.R (only plot the negatives x's) (Use MaxInt Population)
pop <- pop_MaxInt.Gp
x <- pop[pop$Time.Delta.Delay <= 1000 & pop$Death == 1,]$Time.Delta.Delay
y <- pop[pop$Time.Delta.Delay <= 1000 & pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Time.Delta.Delay <= 1000 & pop$Death == 2,]$Time.Delta.Delay
y2 <- pop[pop$Time.Delta.Delay <= 1000 & pop$Death == 2,]$Int.Max.R
# cols <- getDensityColors(x,y)
xcol <- 'black'
ycol <- R
xlim = c(min(c(x,x2)),max(c(x,x2)))
ylim = c(min(c(y,y2)),max(c(y,y2)))
pdf('Figure_11_a.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1,1)
paperPlot(letter='b',xlimit=xlim, ylimit=ylim, xlabel='Host Reaction-Time [h]', ylabel='(Virus) Max-Intensity [au]', xcol=xcol, ycol=ycol)
points(x, y, pch=21, col=rgb(0,0,0,0), bg=rgb(0,0,0,0.3), cex=0.85)
points(x2, y2, pch=21, col=rgb(0,0,0,0), bg=rgb(0,0,0,0.3), cex=0.85)
abline(v=0,lty=4)
# text(mean(xlim), max(ylim), labels='when < 0, host activation\ndetected prior to infection', col = 'darkgreen', cex = 1, adj = c(0.5,1))
# arrows(x1=mean(xlim)-.1*(max(xlim)-min(xlim)),y0=0.75*max(ylim),x0=mean(xlim)+.1*(max(xlim)-min(xlim)), lwd=3, col='darkgreen')
dev.off()

######### Test A ##########
duh <- summary(lm(y~x,data.frame(x=c(x,x2),y=log10(c(y,y2)))))
duh
duh <- cor.test(c(x,x2),log10(c(y,y2)))
duh
duh <- my.prop.test(xThresh=0,x1=x,x2=x2,y1=y,y2=y2)
duh

######### Plot B ##########

# b) x=Time.Delta.Delay y=Int.Max.G (only plot the positive x's)
pop <- pop_MaxInt.Gp
x <- pop[pop$Time.Delta.Delay >= -1000 & pop$Death == 1,]$Time.Delta.Delay
y <- pop[pop$Time.Delta.Delay >= -1000 & pop$Death == 1,]$Int.Max.G
x2 <- pop[pop$Time.Delta.Delay >= -1000 & pop$Death == 2,]$Time.Delta.Delay
y2 <- pop[pop$Time.Delta.Delay >= -1000 & pop$Death == 2,]$Int.Max.G
# cols <- getDensityColors(x,y)
xcol <- 'black'
ycol <- G
xlim = c(min(c(x,x2)),max(c(x,x2)))
ylim = c(min(c(y,y2)),max(c(y,y2)))
pdf('Figure_11_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1,1)
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel='Host Reaction-Time [h]', ylabel='(Host) Max-Intensity [au]', xcol=xcol, ycol=ycol)
points(x, y, pch=21, col=rgb(0,0,0,0), bg=rgb(0,0,0,0.3), cex=0.85)
points(x2, y2, pch=21, col=rgb(0,0,0,0), bg=rgb(0,0,0,0.3), cex=0.85)
abline(v=0,lty=4)
# text(mean(xlim), max(ylim), labels='when > 0, host activation\ndetected after infection', col = 'darkred', cex = 1, adj = c(0.5,1))
# arrows(x0=mean(xlim)-.1*(max(xlim)-min(xlim)),y0=0.75*max(ylim),x1=mean(xlim)+.1*(max(xlim)-min(xlim)), lwd=3, col='darkred')
dev.off()

######### Test B ##########
duh <- summary(lm(y~x,data.frame(x=c(x,x2),y=log10(c(y,y2)))))
duh
duh <- cor.test(x=c(x,x2),y=log10(c(y,y2)))
duh
duh <- my.prop.test(xThresh=0,x1=x,x2=x2,y1=y,y2=y2)
duh

######### Plot C ##########

# c) (Use the Max.Int population here)
pop <- pop_MaxInt.Gp
pop$Ratio.Max <- pop$Int.Max.R/pop$Int.Max.G
pop$Ratio.Max.Log <- log10(pop$Ratio.Max)
xlim <- c(min(pop$Time.Delta.Delay),max(pop$Time.Delta.Delay))
ylim <- c(min(pop$Ratio.Max),max(pop$Ratio.Max))
pdf('Figure_11_c.pdf', height=0.9*H, width=1.8*W)
paperParams(1,1, labsize=1)
par(mar = c(3.7,4.7,0.8,0.8)); # beyond plot frame [b,l,t,r]
par(mgp = c(2.5,0.8,0)); # placement of axis labels [1], numbers [2] and symbols[3]
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel='Host Reaction-Time [h]', ylabel='(Virus Max-Int.) / (Host Max-Int.)', xcol='black', ycol='black', log='y', yticks=c('0.001', '0.01', '0.01', '0.1', '1', '10', '100', '1000'))
points(pop[pop$Death == 1,]$Time.Delta.Delay, pop[pop$Death == 1,]$Ratio.Max, pch=21, cex=0.7, col='black', bg='black')
points(pop[pop$Death == 2,]$Time.Delta.Delay, pop[pop$Death == 2,]$Ratio.Max, pch=21, cex=0.7, col='red', bg='red')
X.1 <- seq(from=-15, to=0, by=0.5)
X.2 <- seq(from=0, to=15, by=0.5)
# dev.off()
library(drc)
fitResults <- data.frame(X_1=numeric(0), Y_1=numeric(0), X_2=numeric(0), Y_2=numeric(0), SSE=numeric(0))
for(x1 in X.1)
{
    for(x2 in X.2)
    {
        defs <- DefinePWFct(x1,x2)
        pNames.pw <- defs$a
        fct.pw <- defs$b
        ssfct.pw <- defs$c
        fctDef.pw <- defs$d
        pwModel <- drm(Ratio.Max.Log~Time.Delta.Delay,data=pop,pmodels=data.frame(1,1),fct=fctDef.pw,na.action=na.omit,separate=FALSE, type='continuous')
        # Gather results from the fitting into a more accesible table format
        temp <- coef(pwModel) # grab the fitted parameters from the model object
        fitResults <- rbind(fitResults, data.frame(X_1 = x1, Y_1 = temp[1], X_2 = x2, Y_2 = temp[2], SSE = sum(residuals(pwModel)^2)))
    }
}

# library(scatterplot3d)
# scatterplot3d(fitResults$X_1, fitResults$X_2, fitResults$SSE, pch=20, highlight.3d=TRUE)#color=rgb(0,0,0,0), bg=rgb(0,0,0,0.5+fitResults$X_2/(2*min(fitResults$X_3))))

bestFit <- which(fitResults$SSE == min(fitResults$SSE))
bestFit <- fitResults[bestFit,]
print(bestFit)

# Reset starter function to be near gloabl min
startDaFunc <- function(data) 
{ 
    return(c(bestFit$Y_1, bestFit$Y_2, bestFit$X_1, bestFit$X_2))
    #     return(c(-4.711, -0.9105, 6, 1.181)) 
}
daFuncDef <- list(name='pwlin',fct=daFunc, ssfct=startDaFunc, names=c('X_1','Y_1','X_2','Y_2'));

# Fit all parameters to refine optimization
pwModel <- drm(Ratio.Max.Log~Time.Delta.Delay,data=pop,pmodels=data.frame(1,1,1,1),fct=daFuncDef,na.action=na.omit,separate=FALSE, type='continuous')
# duh <- daFunc2(pop$Time.Delta.Delay,coef(pwModel))
# duh2 <- daFunc2(pop$Time.Delta.Delay,startDaFunc())
# plot(pop$Time.Delta.Delay, pop$Ratio.Max.Log, type='p', col=rgb(0,0,0,0.5))
# points(pop$Time.Delta.Delay, duh, col='red')
# points(pop$Time.Delta.Delay, duh2, col='green')
summary(pwModel)

# plot fit results and mean indicators over data
from <- -40
to <- 40
lines(cbind(c(from,bestFit$X_1,bestFit$X_2,to),10^(c(bestFit$Y_1,bestFit$Y_1,bestFit$Y_2,bestFit$Y_2))))
lines(cbind(c(from,coef(bestFit)['x1:(Intercept)'],coef(bestFit)['x2:(Intercept)'],to),10^(c(coef(bestFit)['y1:(Intercept)'],coef(bestFit)['y1:(Intercept)'],coef(bestFit)['y2:(Intercept)'],coef(bestFit)['y2:(Intercept)']))))
abline(h=10^(mean(pop$Ratio.Max.Log)),lty=4)
abline(v=mean(pop$Time.Delta.Delay),lty=4)
abline(v=0,lty=1)

text(0-1, 5*min(ylim), labels='host activation\ndetected prior to infection', col = 'darkgreen', cex = 1, adj = c(1,1))
text(0+1, 5*min(ylim), labels='host activation\ndetected after infection', col = 'darkred', cex = 1, adj = c(0,1))
arrows(x0=-1,y0=0.03,x1=-3, lwd=3, col='darkgreen')
arrows(x0=1,y0=0.03, x1=3, lwd=3, col='darkred')

1-bestFit$SSE/sum((pop$Ratio.Max.Log-mean(pop$Ratio.Max.Log))^2)

# linFit <- lm(Ratio.Max.Log~Time.Delta.Delay, data=pop, na.action=na.omit)
# abline(coef(linFit))

# Finish the plot
dev.off()



