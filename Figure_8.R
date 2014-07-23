####################################
########## Figure 2-4     ##########
####################################

# Also compile the results to compare GFP+ and GFP- populations

activationSummary <- data.frame(Mean=numeric(), p=numeric())

########## Red params only #########

########## GFP + ##########
# Max.R vs Alpha.R (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp)
lx <- pop$NLFit.Alpha.R
ly <- pop$Int.Max.R # Don't filter by death because we want the regression to be on red and black points
x <- pop[pop$Death == 1,]$NLFit.Alpha.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.R # Lysed
y2 <- pop[pop$Death == 2,]$Int.Max.R # Lysed
xcol <- R
ycol <- R
xlim = xlim.Alpha
ylim = xlim.Int.Max.R
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_8_a.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='a',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' [1/h]')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.R vs Time.Rise.R (just use maxInt pop)
pop <- pop_MaxInt.Gp
lx <- pop$Time.Rise.R
ly <- pop$Int.Max.R
x <- pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$Int.Max.R
xcol <- R
ycol <- R
xlim = xlim.Time.Rise
ylim = xlim.Int.Max.R
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
Rcor <- cor.test(lx, ly)
pdf('Figure_8_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='b',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.R vs Alpha.R*Time.Rise.R  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp)
lx <- pop$NLFit.Alpha.R*pop$Time.Rise.R
ly <- pop$Int.Max.R
x <- pop[pop$Death == 1,]$NLFit.Alpha.R*pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.R*pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$Int.Max.R
xcol <- R
ycol <- R
xlim = c(min(lx),max(lx))
ylim = xlim.Int.Max.R
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_8_c.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' * Rise-Time')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Alpha.R vs Time.Rise.R  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gp, pop_Alpha.Gp)
lx <- pop$Time.Rise.R
ly <- pop$NLFit.Alpha.R
x <- pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$NLFit.Alpha.R
x2 <- pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$NLFit.Alpha.R
xcol <- R
ycol <- R
xlim = xlim.Time.Rise
ylim = xlim.Alpha
fit <- lm(y~x, data.frame(x=lx, y=ly))
Rcor <- cor.test(lx, ly)
pdf('Figure_8_d.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='d',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel=expression(paste('', alpha, ' [1/h]')), xcol=xcol, ycol=ycol)
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve((coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

########## GFP - ##########
# Max.R vs Alpha.R (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gn, pop_Alpha.Gn)
lx <- pop$NLFit.Alpha.R
ly <- pop$Int.Max.R # Don't filter by death because we want the regression to be on red and black points
x <- pop[pop$Death == 1,]$NLFit.Alpha.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.R # Lysed
y2 <- pop[pop$Death == 2,]$Int.Max.R # Lysed
xcol <- R
ycol <- R
xlim = xlim.Alpha
ylim = xlim.Int.Max.R
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_8_e.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='e',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' [1/h]')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.R vs Time.Rise.R (just use maxInt pop)
pop <- pop_MaxInt.Gn
lx <- pop$Time.Rise.R
ly <- pop$Int.Max.R
x <- pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$Int.Max.R
xcol <- R
ycol <- R
xlim = xlim.Time.Rise
ylim = xlim.Int.Max.R
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_8_f.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='f',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Max.R vs Alpha.R*Time.Rise.R  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gn, pop_Alpha.Gn)
lx <- pop$NLFit.Alpha.R*pop$Time.Rise.R
ly <- pop$Int.Max.R
x <- pop[pop$Death == 1,]$NLFit.Alpha.R*pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.R*pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$Int.Max.R
xcol <- R
ycol <- R
xlim = c(min(lx),max(lx))
ylim = xlim.Int.Max.R + 1
fit <- lm(y~x, data.frame(x=lx, y=log10(ly)))
Rcor <- cor.test(lx, log10(ly))
pdf('Figure_8_g.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='g',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('', alpha, ' * Rise-Time')), ylabel='Max-Intensity [au]', xcol=xcol, ycol=ycol, log='y')
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve(10^(coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()

# Alpha.R vs Time.Rise.R  (Intersect MaxInt and Alpha pops)
pop <- merge(pop_MaxInt.Gn, pop_Alpha.Gn)
lx <- pop$Time.Rise.R
ly <- pop$NLFit.Alpha.R
x <- pop[pop$Death == 1,]$Time.Rise.R
y <- pop[pop$Death == 1,]$NLFit.Alpha.R
x2 <- pop[pop$Death == 2,]$Time.Rise.R
y2 <- pop[pop$Death == 2,]$NLFit.Alpha.R
xcol <- R
ycol <- R
xlim = xlim.Time.Rise
ylim = xlim.Alpha
fit <- lm(y~x, data.frame(x=lx, y=(ly)))
Rcor <- cor.test(lx, (ly))
pdf('Figure_8_h.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='h',xlimit=xlim, ylimit=ylim, xlabel='Rise-Time [h]', ylabel=expression(paste('', alpha, ' [1/h]')), xcol=xcol, ycol=ycol)
points(x, y, pch=21, bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
curve((coef(fit)[1] + x*coef(fit)[2]), xlim=xlim, log='y', lty=2, add=TRUE)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(Rcor$estimate,2))), pos=2)
dev.off()