library(foreign)
path = '/Users/jaywarrick/Desktop/SecretionData/'
AMean <- subset(read.arff('/Users/jaywarrick/Desktop/SecretionData/A Mean - x0_y0.arff'), Measurement=='B')
AMedian <- subset(read.arff('/Users/jaywarrick/Desktop/SecretionData/A Median - x0_y0.arff'), Measurement=='B')
BMean <- read.arff('/Users/jaywarrick/Desktop/SecretionData/B Mean - x0_y0.arff')
BMeanSecretion <- subset(BMean, Measurement=='G')
BMeanBead <- subset(BMean, Measurement=='R')
BMedian <- read.arff('/Users/jaywarrick/Desktop/SecretionData/B Median - x0_y0.arff')
BMedianSecretion <- subset(BMedian, Measurement=='G')
BMedianBead <- subset(BMedian, Measurement=='R')

filename=outputPlotPath;
overallCex = 1;
legendFontCex = 0.8;
axisLabelFontCex = 1.2;
pointCex = 0.5;
lineWidth=2;
jpegH=6;
jpegW=6;



#jpeg(filename,width=jpegW,height=jpegH,units='in',res=300)
CellSignal <- AMean$Value - AMedian$Value
SecretionSignal <- BMeanSecretion$Value - BMedianSecretion$Value
BeadSignal <- BMeanBead$Value - BMedianBead$Value
temp1 <- data.frame(ImRow=AMean$ImRow, ImCol=AMean$ImCol, Point=AMean$Point, Cell=CellSignal, bead=BeadSignal, Secretion=SecretionSignal, Signal=SecretionSignal/BeadSignal)

par(cex=overallCex, cex.lab=axisLabelFontCex, cex.axis=axisLabelFontCex, mar=c(4.5,4.5,2.1,2.1), lwd=lineWidth)
plot(CellSignal, SecretionSignal/BeadSignal, pch=21, cex=pointCex, ylim=c(-2,2), type='p', log='x')
secretionBGHist <- hist(BMedianSecretion$Value, breaks=100, xlab='Secretion Channel Background Intensity [au]', main='')
print(secretionBGHist$mids)
positiveIndicies <- BeadSignal > 10 & BMedianSecretion$Value < 675
temp2 <- data.frame(ImRow=AMean$ImRow[positiveIndicies], ImCol=AMean$ImCol[positiveIndicies], Point=AMean$Point[positiveIndicies], Cell=CellSignal[positiveIndicies], bead=BeadSignal[positiveIndicies], Secretion=SecretionSignal[positiveIndicies], Signal=SecretionSignal[positiveIndicies]/BeadSignal[positiveIndicies])


plot(temp2$Cell, temp2$Signal, pch=21, cex=pointCex, ylim=c(-0.5,3), log='x', type='p', xlab='Nuclear Signal [au]', ylab='Secrection Signal / Bead Signal');

CellNoise_0SD <- mad(temp2$Signal[temp2$Cell < 10])
CellNoise_0Median <- median(temp2$Signal[temp2$Cell < 10])
Thresh <- CellNoise_0Median + 5*CellNoise_0SD
CellThresh <- median(temp2$Cell[temp2$Cell<10]) + 4*mad(temp2$Cell[temp2$Cell<10])

print(Thresh)
PercentSecreting <- sum(temp2$Signal[temp2$Cell >= 10] > Thresh)/length(temp2$Cell >= 10)
PercentFalsePositive <- sum(temp2$Signal[temp2$Cell < 10] > Thresh)/length(temp2$Cell < 10)
print(PercentSecreting)
print(PercentFalsePositive)
abline(h=CellNoise_0Median,col='red')
abline(h=Thresh,col='blue');
abline(v=CellThresh,col='blue')

DaSignal <- log(temp2$Signal[temp2$Cell > 10 & temp2$Signal < 3])
DaWells <- log(temp2$Signal[temp2$Cell < 10 & temp2$Signal < 3])
histSignal <- hist(DaSignal, breaks=30, plot=T)
histWells <- hist(DaWells, breaks=30, plot=T)
plot(histSignal$x, histSignal$mids)



pdf(paste(path,"ForSecretion",'.pdf',sep=''), width=4, height=3)
par(mar=c(4,4,1,1)+0.1)
histX1 <- histWells$mids
histY1 <- histWells$density
plot(histX1, histY1, type='l', xlab='ln(IL8 Signal Intensity per Bead) [au]', ylab='Normalized Frequency', xlim=c(-8,2), ylim=c(0,0.4), add=T)
polygon(pad(histX1, zeros=FALSE), pad(histY1, zeros=TRUE),col=rgb(0,0,0.7,0.7))
histX2 <- histSignal$mids
histY2 <- histSignal$density
lines(histX2, histY2, type='l', xlim=c(-8,2))
polygon(pad(histX2, zeros=FALSE), pad(histY2, zeros=TRUE),col=rgb(0.7,0,0,0.7))
legend('topright', legend=c('Bead Background Levels','MDA-MB231 Secretion'), fill=c(rgb(0,0,0.7,0.7), rgb(0.7,0,0,0.7)), cex=0.5)
dev.off()
