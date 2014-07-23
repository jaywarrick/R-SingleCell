# Alarid / Lang ER Distributions

path <- '/Users/jaywarrick/Documents/MMB_Lab_2/Grants/2014 - PSOC/Elaine/'
cell <- 'HS5'
setwd(path)
cocultureFile <- paste(path, cell, '.txt', sep='');
monoFile <- paste(path, 'Mono.txt', sep='');

coculture <- read.table(cocultureFile, header=TRUE, sep='\t')
mono <- read.table(monoFile, header=TRUE, sep='\t')

metric <- 'Er'

logCoculture <- log(coculture[,metric])
logMono <- log(mono[,metric])
x <- logCoculture
cocultureGuess <- c(mean(x), sd(x))
cocultureDensity <- density(x)

x <- logMono
monoGuess <- c(mean(x), sd(x))
monoGuess2 <- c(3,0.5,4.8,0.3,.7)

f <- function(par,x)
{
    m <- par[1]
    sd <- par[2]
    rhat <- dnorm(x=x, mean=m, sd=sd)
    return(rhat)
}

f2 <- function(par,x)
{
    m1 <- par[1]
    sd1 <- par[2]
    m2 <- par[3]
    sd2 <- par[4]
    r <- par[5]
    rhat <- r * dnorm(x=x, mean=m1, sd=sd1) + (1-r) * dnorm(x=x, mean=m2, sd=sd2)
    return(rhat)
}

fsse <- function(par,data)
{
    d <- density(data)
    r <- d$y
    x <- d$x
    rhat <- f(par,x)
    sum((r - rhat)^2)
}

fsse2 <- function(par,data)
{
    d <- density(data)
    r <- d$y
    x <- d$x
    rhat <- f2(par,x)
    sum((r - rhat)^2)
}

pad <- function(x, zeros=TRUE)
{
    if(zeros)
    {
        return(c(0,x,0,0))
    }
    else
    {
        return(c(min(x),x,max(x),min(x)))
    }
}

y <- logCoculture
cocultureFit <- optim(cocultureGuess, fsse, method="BFGS", data = y, control=list(reltol=1e-10))
yHist <- hist(y, breaks=25, plot=FALSE)
plot(yHist$mids, yHist$density, type='l')
histX <- yHist$mids
histY <- yHist$density
polygon(pad(histX, zeros=FALSE), pad(histY, zeros=TRUE),col=gray(0.7))
lines(density(y)$x, f(cocultureFit$par, density(y)$x), col='black', lwd=2)
cocultureFit

# Test reverse order
monoFit <- optim(monoGuess, fsse, method="BFGS", data = mono[,metric], control=list(reltol=1e-10))
yHist2 <- hist(mono[,metric], breaks=40, plot=FALSE)
plot(log(yHist2$mids), yHist2$density, type='l')
histX2 <- log(yHist2$mids)
histY2 <- yHist2$density
polygon(pad(histX2, zeros=FALSE), pad(histY2, zeros=TRUE),col=gray(0.7))
lines(log(density(mon[,metric])$x, f(monoGuess, density(mono[,metric])$x), col='black', lwd=2)
monoFit

y2 <- logMono
monoFit <- optim(monoGuess, fsse, method="BFGS", data = y2, control=list(reltol=1e-10))
yHist2 <- hist(y2, breaks=25, plot=FALSE)
plot(yHist2$mids, yHist2$density, type='l')
histX2 <- yHist2$mids
histY2 <- yHist2$density
polygon(pad(histX2, zeros=FALSE), pad(histY2, zeros=TRUE),col=gray(0.7))
lines(density(y2)$x, f(monoGuess, density(y2)$x), col='black', lwd=2)
monoFit
monoFit2 <- optim(monoGuess2, fsse2, method="BFGS", data = y2, control=list(reltol=1e-10))
yHist2b <- hist(y2, breaks=25, plot=FALSE)
# pdf(paste(path,"ForElaine",'.pdf',sep=''), width=4, height=3)
par(mar=c(4,4,1,1)+0.1)
plot(yHist2b$mids, yHist2b$density, type='l', lwd=0.5, ylab='Normalized Frequency', xlab='ln(ER Integrated Intensity)')
histX2b <- yHist2b$mids
histY2b <- yHist2b$density
polygon(pad(histX2b, zeros=FALSE), pad(histY2b, zeros=TRUE),col=rgb(0,0,0.7,0.7))
lines(density(y2)$x, f2(monoFit2$par, density(y2)$x), col='black', lwd=2)
lines(density(y2)$x, monoFit2$par[5]*f(monoFit2$par[1:2], density(y2)$x), col='black', lwd=2, lty=2)
lines(density(y2)$x, (1-monoFit2$par[5])*f(monoFit2$par[3:4], density(y2)$x), col='black', lwd=2, lty=2)
monoFit2

# Energy Calculation
ER <- density(y2)$x
energy <- -0.5*log(f2(monoFit2$par, ER))
# pdf(paste(path,"EnergyForElaine",'.pdf',sep=''), width=4, height=3)
par(mar=c(3,4,1,1)+0.1)
plot(ER, energy, col='black', type='l',xlab='', ylab='', lwd=2)
# plot(yHist2b$mids, yHist2b$density, type='l', lwd=0.5, ylab='Normalized Frequency', xlab='ln(ER Integrated Intensity)')
histX2b <- yHist2b$mids
histY2b <- yHist2b$density
polygon(pad(histX2b, zeros=FALSE), pad(histY2b, zeros=TRUE),col=rgb(0,0,0.7,0.7))
lines(ER, energy, col='black', lwd=2)
mtext(expression(paste(Phi,'/D',sep='')), side=2, line=2)
mtext('ER Intensity [au]', side=3, line=2)
# dev.off()

daMax = max(energy[ER > 3 & ER < 4])
lMin = min(energy)
rMin = min(energy[ER > 4.2])
dMax
lMin
rMin

ER

pdf(paste(path,cell,'.pdf',sep=''), width=4, height=3)
xRange = range(histX, histX2)
yRange = range(histY, histY2)
par(mar=c(4,4,1,1)+0.1)
plot(yHist$mids, yHist$density, type='l', xlim=xRange, ylim=yRange, ylab='Normalized Frequency', xlab='ln(ER Integrated Intensity)', lty=0)
histX <- yHist$mids
histY <- yHist$density
polygon(pad(histX, zeros=FALSE), pad(histY, zeros=TRUE),col=rgb(0.7,0,0,0.5), border=FALSE)
lines(density(y)$x, f(cocultureFit$par, density(y)$x), col='black', lwd=2)
lines(yHist2$mids, yHist2$density, lty=0)
histX2 <- yHist2$mids
histY2 <- yHist2$density
polygon(pad(histX2, zeros=FALSE), pad(histY2, zeros=TRUE),col=rgb(0,0,0.7,0.5), border=FALSE)
# lines(density(y2)$x, f(monoFit$par, density(y2)$x), col='black', lwd=2, lty=2)
# dev.off()

y2 <- logMono
# monoFit <- optim(monoGuess, fsse, method="BFGS", data = y2, control=list(reltol=1e-10))
# yHist2 <- hist(y2, breaks=25, plot=FALSE)
# plot(yHist2$mids, yHist2$density, type='l')
# histX2 <- yHist2$mids
# histY2 <- yHist2$density
# polygon(pad(histX2, zeros=FALSE), pad(histY2, zeros=TRUE),col=gray(0.7))
# lines(density(y2)$x, f(monoGuess, density(y2)$x), col='black', lwd=2)
# monoFit
monoFit2 <- optim(monoGuess2, fsse2, method="BFGS", data = y2, control=list(reltol=1e-10))
yHist2b <- hist(y2, breaks=25, plot=FALSE)
# pdf(paste(path,"ForElaine",'.pdf',sep=''), width=4, height=3)
par(mar=c(4,4,1,1)+0.1)
# lines(yHist2b$mids, yHist2b$density, type='l', lwd=0.5, ylab='Normalized Frequency', xlab='ln(ER Integrated Intensity)')
histX2b <- yHist2b$mids
histY2b <- yHist2b$density
# polygon(pad(histX2b, zeros=FALSE), pad(histY2b, zeros=TRUE),col=rgb(0,0,0.7,0.7))
lines(density(y2)$x, f2(monoFit2$par, density(y2)$x), col='black', lty=2, lwd=2)
# lines(density(y2)$x, monoFit2$par[5]*f(monoFit2$par[1:2], density(y2)$x), col='black', lwd=2, lty=2)
# lines(density(y2)$x, (1-monoFit2$par[5])*f(monoFit2$par[3:4], density(y2)$x), col='black', lwd=2, lty=2)
dev.off()
monoFit2



