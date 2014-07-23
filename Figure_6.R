####################################
######### Figure 2-3       #########
####################################

# Max.Int.R
pop <- pop_MaxInt.All
y1 <- pop$Int.Max.R
y2 <- pop$Int.Max.G
lims <- range(y1[y1 > 0], y2[y2 > 0], na.rm=TRUE)
breaks <- seq(from=0, to=lims[2], length.out=10)
hist1 <- hist(y1, breaks=breaks, plot=FALSE)
hist5 <- hist(y2, breaks=breaks, plot=FALSE)
med1 <- median(y1, na.rm=TRUE)
mad1 <- mad(y1, na.rm=TRUE)
med5 <- median(y2, na.rm=TRUE)
mad5 <- mad(y2, na.rm=TRUE)
range1 <- range(y1[y1 > 0], na.rm=TRUE)
range5 <- range(y2[y2 > 0], na.rm=TRUE)
mu1 <- mean(y1, na.rm=TRUE)
mu5 <- mean(y2, na.rm=TRUE)
xlim.Int.Max.R <- lims
xlim.Int.Max.G <- lims
ylim.Int.Max <- c(0,max(hist1$density,hist5$density))

# Time.Delay.R
pop <- pop_Delay.All
y1 <- pop$Time.Delay.R
y2 <- pop$Time.Delay.G
lims <- range(y1, y2, na.rm=TRUE)
breaks <- seq(from=lims[1], to=lims[2], length.out=10)
hist2 <- hist(y1, breaks=breaks, plot=FALSE)
hist6 <- hist(y2, breaks=breaks, plot=FALSE)
med2 <- median(y1, na.rm=TRUE)
mad2 <- mad(y1, na.rm=TRUE)
med6 <- median(y2, na.rm=TRUE)
mad6 <- mad(y2, na.rm=TRUE)
range2 <- range(y1, na.rm=TRUE)
range6 <- range(y2, na.rm=TRUE)
mu2 <- mean(y1, na.rm=TRUE)
mu6 <- mean(y2, na.rm=TRUE)
xlim.Time.Delay <- lims
ylim.Time.Delay <- c(0,max(hist2$density,hist6$density))

# Max.Alpha.R
pop <- pop_MaxInt.All
y1 <- pop$NLFit.Alpha.R
y2 <- pop$NLFit.Alpha.G
lims <- range(y1, y2, na.rm=TRUE)
breaks <- seq(from=lims[1], to=lims[2], length.out=10)
hist3 <- hist(y1, breaks=breaks, plot=FALSE)
hist7 <- hist(y2, breaks=breaks, plot=FALSE)
med3 <- median(y1, na.rm=TRUE)
mad3 <- mad(y1, na.rm=TRUE)
med7 <- median(y2, na.rm=TRUE)
mad7 <- mad(y2, na.rm=TRUE)
range3 <- range(y1, na.rm=TRUE)
range7 <- range(y2, na.rm=TRUE)
mu3 <- mean(y1, na.rm=TRUE)
mu7 <- mean(y2, na.rm=TRUE)
xlim.Alpha <- lims
ylim.Alpha <- c(0,max(hist3$density,hist7$density))

# Time.Rise
pop <- pop_MaxInt.All
y1 <- pop$Time.Rise.R
y2 <- pop$Time.Rise.G
lims <- range(y1, y2, na.rm=TRUE)
breaks <- seq(from=lims[1], to=lims[2], length.out=10)
hist4 <- hist(y1, breaks=breaks, plot=FALSE)
hist8 <- hist(y2, breaks=breaks, plot=FALSE)
med4 <- median(y1, na.rm=TRUE)
mad4 <- mad(y1, na.rm=TRUE)
med8 <- median(y2, na.rm=TRUE)
mad8 <- mad(y2, na.rm=TRUE)
range4 <- range(y1, na.rm=TRUE)
range8 <- range(y2, na.rm=TRUE)
mu4 <- mean(y1, na.rm=TRUE)
mu8 <- mean(y2, na.rm=TRUE)
xlim.Time.Rise <- lims
ylim.Time.Rise <- c(0,max(hist4$density,hist8$density))

# rm(histogram)
# rm(xlabel)
# rm(xlim)
# rm(ylim)
# rm(ylabel)
# rm(xcol)
# rm(bar.col)
# histogram <- hist4
# xlabel <- 'Rise-Time [h]'
# xlim <- xlim.Time.Rise
# ylim <- ylim.Time.Rise
# ylabel <- ''
# xcol <- R
# bar.col <- R

# composite plot
pdf('Figure_6_a.pdf', height=0.9*H, width=0.9*W) # Height used to be 0.6*H
paperParams(1,1)
#paperHist2('a',hist3, hist7, xlabel='Production-Rate [1/h]', xlim=xlim.Alpha, ylim=ylim.Alpha)
paperHist('a',hist3, hist7, xlabel='Production-Rate [1/h]', xlim=xlim.Alpha, ylim=ylim.Alpha, ylabel='', xcol=R, bar.col=R)
med <- med3
mad <- mad3
range <- range3
mu <- mu3
n <- length(hist3$mids)
text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
dev.off()

pdf('Figure_6_b.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
#paperHist2('b',hist2, hist6, xlabel='Delay-Time [h]', xlim=xlim.Time.Delay, ylim=ylim.Time.Delay, ylabel='')
paperHist('b',hist2, xlabel='Delay-Time [h]', xlim=xlim.Time.Delay, ylim=ylim.Time.Delay, ylabel='', xcol=R, bar.col=R)
med <- med2
mad <- mad2
range <- range2
mu <- mu2
text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
dev.off() # turn off figure device

pdf('Figure_6_c.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
#paperHist2('c',hist4, hist8, xlabel='Rise-Time [h]', xlim=xlim.Time.Rise, ylim=ylim.Time.Rise, ylabel='')
paperHist('c',hist4, xlabel='Rise-Time [h]', xlim=xlim.Time.Rise, ylim=ylim.Time.Rise, ylabel='', xcol=R, bar.col=R)
med <- med4
mad <- mad4
range <- range4
mu <- mu4
text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
dev.off() # turn off figure device

pdf('Figure_6_d.pdf', height=0.9*H, width=0.9*W)
paperParams(1,1)
#paperHist2('d',hist1, hist5, xlabel='Max-Intensity [au]', xlim=xlim.Int.Max.R, ylim=ylim.Int.Max, ylabel='')
#legend('right',legend=c('Virus','Host'), fill=c(rgb(1,0,0,0.5),rgb(0,1,0,0.2)), yjust=0)
paperHist('d',hist1, xlabel='Max-Intensity [au]', xlim=xlim.Int.Max.R, ylim=ylim.Int.Max, ylabel='', xcol=R, bar.col=R)
med <- med1
mad <- mad1
range <- range1
mu <- mu1

text(max(xlim.Int.Max.R),max(ylim.Int.Max),adj=c(1,1),labels=bquote(bar(x): .(round(mu,0))), cex=0.7)
text(max(xlim.Int.Max.R),max(ylim.Int.Max),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,0))), cex=0.7)
text(max(xlim.Int.Max.R),max(ylim.Int.Max),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,0))), cex=0.7)
text(max(xlim.Int.Max.R),max(ylim.Int.Max),adj=c(1,4.7),labels=bquote(range: .(round(range[1],0))-.(round(range[2],0))), cex=0.7)
dev.off() # turn off figure device

# pdf('Figure_6_e.pdf', height=0.6*H, width=0.9*W)
# paperParams(1,1)
# paperHist('e',hist7, xlabel='Production-Rate [1/h]', xlim=xlim.Alpha, ylim=ylim.Alpha, xcol=G, bar.col=G)
# med <- med7
# mad <- mad7
# range <- range7
# mu <- mu7
# text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
# text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
# text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
# text(max(xlim.Alpha),max(ylim.Alpha),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
# dev.off() # turn off figure device
# 
# pdf('Figure_6_f.pdf', height=0.6*H, width=0.9*W)
# paperParams(1,1)
# paperHist('f',hist6, xlabel='Delay-Time [h]', xlim=xlim.Time.Delay, ylim=ylim.Time.Delay, ylabel='', xcol=G, bar.col=G)
# med <- med6
# mad <- mad6
# range <- range6
# mu <- mu6
# text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
# text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
# text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
# text(max(xlim.Time.Delay),max(ylim.Time.Delay),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
# dev.off() # turn off figure device
# 
# pdf('Figure_6_g.pdf', height=0.6*H, width=0.9*W)
# paperParams(1,1)
# paperHist('g',hist8, xlabel='Rise-Time [h]', xlim=xlim.Time.Rise, ylim=ylim.Time.Rise, ylabel='', xcol=G, bar.col=G)
# med <- med8
# mad <- mad8
# range <- range8
# mu <- mu8
# text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,1),labels=bquote(bar(x): .(round(mu,1))), cex=0.7)
# text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,1))), cex=0.7)
# text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,1))), cex=0.7)
# text(max(xlim.Time.Rise),max(ylim.Time.Rise),adj=c(1,4.7),labels=bquote(range: .(round(range[1],1))-.(round(range[2],1))), cex=0.7)
# dev.off() # turn off figure device
# 
# pdf('Figure_6_h.pdf', height=0.6*H, width=0.9*W)
# paperParams(1,1)
# paperHist('h',hist5, xlabel='Max-Intensity [au]', xlim=xlim.Int.Max.G, ylim=ylim.Int.Max, ylabel='', xcol=G, bar.col=G)
# med <- med5
# mad <- mad5
# range <- range5
# mu <- mu5
# text(max(xlim.Int.Max.G),max(ylim.Int.Max),adj=c(1,1),labels=bquote(bar(x): .(round(mu,0))), cex=0.7)
# text(max(xlim.Int.Max.G),max(ylim.Int.Max),adj=c(1,2.5),labels=bquote(tilde(x): .(round(med,0))), cex=0.7)
# text(max(xlim.Int.Max.G),max(ylim.Int.Max),adj=c(1,4.2),labels=bquote(sigma: .(round(mad,0))), cex=0.7)
# text(max(xlim.Int.Max.G),max(ylim.Int.Max),adj=c(1,4.7),labels=bquote(range: .(round(range[1],0))-.(round(range[2],0))), cex=0.7)
# dev.off() # turn off figure device
