rm(list=ls())
library(plyr)
library(zoo)
library(foreign)
library(tiff)
library(data.table)
source('/Users/jaywarrick/Public/Dropbox/GitHub/R-SingleCell/preprocessingHelperFunctions.R')

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

# preprocess and write the preprocessed data to files
temp <- preprocessVirusData(data, virusType='M51R')
write.table(temp$single, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Data_1Cell.txt', row.names=FALSE)
write.table(temp$zero, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Data_0Cell.txt', row.names=FALSE)
write.table(temp$thresholds, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Thresholds.txt', row.names=FALSE)
temp <- preprocessVirusData(data, virusType='N1')
write.table(temp$single, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Data_1Cell.txt', row.names=FALSE)
write.table(temp$zero, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Data_0Cell.txt', row.names=FALSE)
write.table(temp$thresholds, '/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/N1_Thresholds.txt', row.names=FALSE)

