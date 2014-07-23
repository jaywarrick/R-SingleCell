##############################
###### Figure 2-5   ##########
##############################

########## Green params only #########

########## GFP + ##########
# Max.G vs Alpha.G (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp.G)
lx <- pop$NLFit.Alpha.G
ly <- pop$Int.Max.G # Don't filter by death because we want the regression to be on red and black points
x <- pop[pop$Death == 1,]$NLFit.Alpha.G
y <- pop[pop$Death == 1,]$Int.Max.G
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.G # Lysed
y2 <- pop[pop$Death == 2,]$Int.Max.G # Lysed
xcol <- G
ycol <- G
xlim = xlim.Alpha
ylim = xlim.Int.Max.G
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_9_a.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='a',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' [1/h]')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.G vs Time.Rise.G (just use maxInt pop)
pop <- pop_MaxInt.Gp
lx <- pop$Time.Rise.G
ly <- pop$Int.Max.G
x <- pop[pop$Death == 1,]$Time.Rise.G
y <- pop[pop$Death == 1,]$Int.Max.G
x2 <- pop[pop$Death == 2,]$Time.Rise.G
y2 <- pop[pop$Death == 2,]$Int.Max.G
xcol <- G
ycol <- G
xlim = xlim.Time.Rise
ylim = xlim.Int.Max.G
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_9_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='b',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.G vs Alpha.G*Time.Rise.G  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp.G)
lx <- pop$NLFit.Alpha.G*pop$Time.Rise.G
ly <- pop$Int.Max.G
x <- pop[pop$Death == 1,]$NLFit.Alpha.G*pop[pop$Death == 1,]$Time.Rise.G
y <- pop[pop$Death == 1,]$Int.Max.G
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.G*pop[pop$Death == 2,]$Time.Rise.G
y2 <- pop[pop$Death == 2,]$Int.Max.G
xcol <- G
ycol <- G
xlim = c(min(lx),max(lx))
ylim = xlim.Int.Max.G
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_9_c.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' * Rise-Time')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Alpha.G vs Time.Rise.G  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp.G)
lx <- pop$Time.Rise.G
ly <- pop$NLFit.Alpha.G
x <- pop[pop$Death == 1,]$Time.Rise.G
y <- pop[pop$Death == 1,]$NLFit.Alpha.G
x2 <- pop[pop$Death == 2,]$Time.Rise.G
y2 <- pop[pop$Death == 2,]$NLFit.Alpha.G
xcol <- G
ycol <- G
fit <- lm(y~x, data.frame(x=lx, y=(ly)))
Rcor <- cor.test(lx, (ly))
xlim = xlim.Time.Rise
ylim = xlim.Alpha
pdf('Figure_9_d.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='d',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel=expression(paste('', alpha, ' [1/h]')), xcol=xcol, ycol=ycol)
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve((coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()


# x <- 1:10
# y <- 10^(0.1*x)
# fit <- lm(y~x, data.frame(x=x, y=log10(y)))
# plot(x,y,log='y')
# cor.test(x, log10(y))
# curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=range(x), log='y', lty=2, add=TRUE)