rm(list=ls())
library(plyr)
library(zoo)
library(foreign)
library(tiff)
library(data.table)
source('/users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/preprocessingHelperFunctions.R')

# Generate Master data files based on arff files from JEX
# masterData <- getData()
# errors <- getErrors()
# data <- masterData$data # Note that data is a data.table not a data.frame
# countData <- masterData$countData # Note that countData is a data.table not a data.frame
# writeData(data, countData, errors)
# andreaData <- getAndreaData()

# Load saved Master data
savedData <- readData()
data <- savedData$data
countData <- savedData$countData
errors <- savedData$errors
rm(savedData)

# Remove ROIS with errors (bubbles etc.)
data <- data[!(data$ID %in% errors$ID),]

# Focus on M51R Data
data <- subset(data, Virus=='M51R')

# Read in the times at which each frame was acquired and fill in the time column with the image frame number (starting at 1)
timedata <- read.csv('/Users/jaywarrick/Desktop/Octave/timeData.dat')
data[, time:=timedata$Time[1:length(Time)], by=ID]

# Calculate background corrected values for each color
data$R.BC <- data$R0_CellMax - data$R0_Mode
data$G.BC <- data$G0_CellMax - data$G0_Mode
data$B.BC <- data$B0_CellMax - data$B0_Mode

# Isolate data with 0 cell in the wells
null <- subset(data, Cell.Count==0)

# Remove wells with more than 1 cell
data <- subset(data, !is.na(Cell.Count))

# Get the mean null signals over time for each ROI and take the median for each
nullSummary <- null[, list(R.BC.Mean=mean(R.BC),G.BC.Mean=mean(G.BC),B.BC.Mean=mean(B.BC)), by=list(Device,Image,ID)] # Time-averaged mean of R, G, and B signals for each ROI with 0-cells in them
nullSummary <- nullSummary[, list(R.BC.Mean.Median=median(R.BC.Mean), G.BC.Mean.Median=median(G.BC.Mean), B.BC.Mean.Median=median(B.BC.Mean)), by=list(Device,Image)] # Median of 0-cell well R, G, and B signals for each image and device

# Subtract the median signal of null wells obtained for image from data in the corresponding images
nullSummary <- data.frame(nullSummary)
data[, R.NULL:=nullSummary[nullSummary$Image==Image[1] & nullSummary$Device==Device[1],]$R.BC.Mean.Median, by=list(Device,Image)]
data[, G.NULL:=nullSummary[nullSummary$Image==Image[1] & nullSummary$Device==Device[1],]$G.BC.Mean.Median, by=list(Device,Image)]
data[, B.NULL:=nullSummary[nullSummary$Image==Image[1] & nullSummary$Device==Device[1],]$B.BC.Mean.Median, by=list(Device,Image)]
data[, R.BC.NULL:=R.BC-R.NULL]
data[, G.BC.NULL:=G.BC-G.NULL]
data[, B.BC.NULL:=B.BC-B.NULL]

# Apply illumination correction
countData <- data.frame(countData)
data[, IF.Factor:=countData[countData$ID==ID[1],]$IF.Factor, by=ID]
data[, R.Final:=R.BC.NULL/IF.Factor]
data[, G.Final:=G.BC.NULL/IF.Factor]
data[, B.Final:=B.BC.NULL/IF.Factor]

# Separate the data into data for 0-cell wells and 1-cell wells
single <- subset(data, Cell.Count==1)
zero <- subset(data, Cell.Count==0)

# Create a quick plot for a sense check
plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$R.Final)))
single[, lines(Time,R.Final, col='red'), by=ID]
zero[, lines(Time,R.Final, col='black'), by=ID]

plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$G.Final)))
single[, lines(Time,G.Final, col='green'), by=ID]
zero[, lines(Time,G.Final, col='black'), by=ID]

plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$B.Final)))
single[, lines(Time,B.Final, col='blue'), by=ID]
zero[, lines(Time,B.Final, col='black'), by=ID]

# Determine the threshold based on 0-cell data
zeroSummary <- zero[, list(R.Mean=mean(R.Final),G.Mean=mean(G.Final),B.Mean=mean(B.Final)), by=list(Device,Image,ID)]
zeroSummary <- zeroSummary[, list(R.Mean.Median=median(R.Mean), G.Mean.Median=median(G.Mean), B.Mean.Median=median(B.Mean), R.Mean.StdDev=sd(R.Mean), G.Mean.StdDev=sd(G.Mean), B.Mean.StdDev=sd(B.Mean)), by=list(Device,Image)]
zeroSummary[, R.Thresh:=R.Mean.Median+3*R.Mean.StdDev]
zeroSummary[, G.Thresh:=G.Mean.Median+3*G.Mean.StdDev]
zeroSummary[, B.Thresh:=B.Mean.Median+3*B.Mean.StdDev]
thresholds <- zeroSummary[, list(R=max(R.Thresh), G=max(G.Thresh), B=max(B.Thresh))]

# Zero points below threshold and spurrious points above threshold (i.e. less than 4 points in a row above 0)
single[,R:=rollmin(R.Final, thresh=thresholds$R),by=ID]
single[,G:=rollmin(G.Final, thresh=thresholds$G),by=ID]
single[,B:=rollmin(B.Final, thresh=thresholds$B),by=ID]

# For each color, mark trajectories that have a signal
hasSignal <- function(piece){if(max(piece)>0){return(TRUE)}else{return(FALSE)}}
single[,Flag.R:=hasSignal(R),by=ID]
single[,Flag.G:=hasSignal(G),by=ID]
single[,Flag.B:=hasSignal(B),by=ID]

# Create a quick plot for a sense check
plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$R)), log='y')
single[, lines(Time,R, col=rgb(1,0,0,0.2)), by=ID]

plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$G)), log='y')
single[, lines(Time,G, col=rgb(0,1,0,0.2)), by=ID]

plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$B)), log='y')
single[, lines(Time,B, col=rgb(0,0,1,0.2)), by=ID]

# Write the preprocessed data to a file
write.table(single, '/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/M51R_Data_1Cell.txt', row.names=FALSE)
write.table(zero, '/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/M51R_Data_0Cell.txt', row.names=FALSE)
write.table(thresholds, '/Users/jaywarrick/GoogleDrive/SingleCell/AndreaMatlabFiles/JayRFiles/M51R_Thresholds.txt', row.names=FALSE)
