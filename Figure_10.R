# Figure 10: (Include Red lysed dots)
# a) & b) histograms if necessary
# c) compare delay times
# d) piece-wise fit of time correlation between ratio of red and green to relative delay (Use MaxInt population)
# e) x=Time.Delay.Delta y=Max.Int.R (only plot the negatives x's)
# f) x=Time.Delay.Delta y=Max.Int.G (only plot the positive x's)
########################

# a)
pop <- pop_Delay.Gp
x <- pop$Time.Delay.G
y <- pop$Time.Delay.R
# x2 <- pop[pop$Death == 2,]$Time.Delay.G
# y2 <- pop[pop$Death == 2,]$Time.Delay.R
cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = xlim.Time.Delay
ylim = xlim.Time.Delay
pdf('Figure_10_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='a',xlimit=xlim, ylimit=ylim, xlabel='(Host) Delay-Time [hpi]', ylabel='(Virus) Delay-Time [hpi]', xcol=xcol, ycol=ycol)
points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
# points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(a=0,b=1, lty=2)
legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

# Time.Rise.G vs Time.Rise.R  (merge(pop_MaxInt.Gp, pop_MaxInt.Gp))
pop <- merge(pop_MaxInt.Gp, pop_MaxInt.Gp)
x <- pop$Time.Rise.G
y <- pop$Time.Rise.R
# x2 <- pop[pop$Death == 2,]$Time.Rise.G
# y2 <- pop[pop$Death == 2,]$Time.Rise.R
cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = range(c(x,y))
ylim = range(c(x,y))
pdf('Figure_10_c.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='c',xlimit=xlim, ylimit=ylim, xlabel='(Host) Rise-Time [h]', ylabel='(Virus) Rise-Time [h]', xcol=xcol, ycol=ycol)
points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
# points(x, y, pch=21, col='black', bg='black', cex=0.7)
# points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(a=0,b=1, lty=2)
legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()

# Time.Peak.G vs Time.Peak.R  (merge(pop_MaxInt.Gp, pop_MaxInt.Gp))
pop <- merge(pop_MaxInt.Gp, pop_MaxInt.Gp)
x <- pop$Time.Peak.G
y <- pop$Time.Peak.R
# x2 <- pop[pop$Death == 2,]$Time.Rise.G
# y2 <- pop[pop$Death == 2,]$Time.Rise.R
cols <- getDensityColors(x,y)
xcol <- G
ycol <- R
xlim = range(c(x,y))
ylim = range(c(x,y))
pdf('Figure_10_a.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
paperPlot(letter='a',xlimit=xlim, ylimit=ylim, xlabel='(Host) Duration [hpi]', ylabel='(Virus) Duration [hpi]', xcol=xcol, ycol=ycol)
points(cols$data$x, cols$data$y, pch=21, col=cols$data$c, bg=cols$data$c, cex=0.7)
# points(x, y, pch=21, col='black', bg='black', cex=0.7)
# points(x2, y2, pch=21, col='red', bg='red', cex=0.7)
abline(a=0,b=1, lty=2)
legend('topleft', legend=as.character(seq(from=min(cols$data$n), to=max(cols$data$n), by=1)), bg=rgb(1,1,1), pch=21, col=cols$colors, pt.bg=cols$colors, pt.cex=0.7, cex=0.7, y.intersp=0.9)
box(lwd=2)
dev.off()
