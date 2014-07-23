##############################
###### Figure 2-7   ##########
##############################

########## Green params only #########

# Int.Max.G vs Int.Max.R  (merge(pop_MaxInt.Gp, pop_MaxInt.Gp))
pop <- merge(pop_MaxInt.Gp, pop_MaxInt.Gp)
x <- pop[pop$Death == 1,]$Int.Max.G
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$Int.Max.G
y2 <- pop[pop$Death == 2,]$Int.Max.R
# cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = xlim.Int.Max.G
ylim = xlim.Int.Max.R
pdf('Figure_12_a.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='a',xlimit=xlim, ylimit=ylim, log='xy', xlabel='(Host) Max-Intensity [au]', ylabel='(Virus) Max-Intensity [au]', xcol=xcol, ycol=ycol)
# points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
points(x, y, pch=21, col='black', bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(lm(y~x, data.frame(x=log10(c(x,x2)),y=log10(c(y,y2)))), lty=2)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(cor.test(log10(c(x,x2)),log10(c(y,y2)))$estimate,2))), pos=2)
# legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

# Alpha.G vs Int.Max.R  (merge(pop_Alpha.Gp.G, pop_MaxInt.Gp))
pop <- merge(pop_Alpha.Gp.G, pop_MaxInt.Gp)
x <- pop[pop$Death == 1,]$NLFit.Alpha.G
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.G
y2 <- pop[pop$Death == 2,]$Int.Max.R
# cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = xlim.Alpha
ylim = xlim.Int.Max.R
pdf('Figure_12_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='b',xlimit=xlim, ylimit=ylim, log='y', xlabel=expression(paste('(Host) ', alpha, ' [1/h]')), ylabel='(Virus) Max-Intensity [au]', xcol=xcol, ycol=ycol)
# points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
points(x, y, pch=21, col='black', bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(lm(y~x, data.frame(x=c(x,x2),y=log10(c(y,y2)))), lty=2)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(cor.test(c(x,x2),log10(c(y,y2)))$estimate,2))), pos=2)
# abline(a=0,b=1, lty=2)
# legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

# Alpha.G vs Time.Rise.R  (merge(pop_Alpha.Gp.G, pop_MaxInt.Gp))
pop <- merge(pop_Alpha.Gp.G, pop_MaxInt.Gp)
x <- pop[pop$Death == 1,]$NLFit.Alpha.G
y <- pop[pop$Death == 1,]$Time.Rise.R
x2 <- pop[pop$Death == 2,]$NLFit.Alpha.G
y2 <- pop[pop$Death == 2,]$Time.Rise.R
# cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = xlim.Alpha
ylim = xlim.Time.Rise
pdf('Figure_12_c.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel=expression(paste('(Host) ', alpha, ' [1/h]')), ylabel='(Virus) Rise-Time [h]', xcol=xcol, ycol=ycol)
# points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
points(x, y, pch=21, col='black', bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(lm(y~x, data.frame(x=c(x,x2),y=c(y,y2))), lty=2)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(cor.test(c(x,x2),c(y,y2))$estimate,2))), pos=2)
# abline(a=0,b=1, lty=2)
# legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

# Time.Rise.G vs Int.Max.R  (merge(pop_MaxInt.Gp, pop_MaxInt.Gp))
pop <- merge(pop_MaxInt.Gp, pop_MaxInt.Gp)
x <- pop[pop$Death == 1,]$Time.Rise.G
y <- pop[pop$Death == 1,]$Int.Max.R
x2 <- pop[pop$Death == 2,]$Time.Rise.G
y2 <- pop[pop$Death == 2,]$Int.Max.R
# cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = xlim.Time.Rise
ylim = xlim.Int.Max.R
pdf('Figure_12_d.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='d',xlimit=xlim, ylimit=ylim, log='y', xlabel='(Host) Rise-Time [h]', ylabel='(Virus) Max-Intensity [au]', xcol=xcol, ycol=ycol)
# points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
points(x, y, pch=21, col='black', bg='black', cex=0.7)
points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(lm(y~x, data.frame(x=c(x,x2),y=log10(c(y,y2)))), lty=2)
text(max(xlim),min(ylim), labels=c(paste('R = ', round(cor.test(c(x,x2),log10(c(y,y2)))$estimate,2))), pos=2)
# abline(a=0,b=1, lty=2)
# legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

