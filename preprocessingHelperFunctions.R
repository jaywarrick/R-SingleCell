getDevice <- function(x, y)
{
    device <- (y * 3 + x)

    return (device)
}

getImage <- function(ImRow, ImCol)
{
    image <- ImCol*4 + (ImRow + 1)

    return (image)
}

getID <- function(device, image, ROI)
{
    deviceImage <- (device-1)*12 + image
    ID <- 1000*deviceImage + ROI

    return (ID)
}

reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE)
{
    library(plyr)
    idCols <- names(data)
    idCols <- idCols[-which(idCols %in% c('Measurement','Value'))]
    newData <- data.frame()
    measurements <- unique(data$Measurement)
    for(m in measurements)
    {
        if(is.na(baseName))
        {
            newColName <- m
            newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
        }
        else
        {
            newColName <- paste(baseName,'.',m, sep='')
            newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
        }

        temp <- subset(data, Measurement==m)
        temp2 <- temp[,idCols]
        temp2[,newColName] <- temp$Value
        if(nrow(newData) == 0)
        {
            newData <- temp2
        }
        else
        {
            newData <- merge(newData, temp2, by=idCols)
        }
    }

    if(convertToNumeric)
    {
        for(n in idCols)
        {
            newData[,n] <- as.numeric(as.character(newData[,n]))
        }
    }

    return(newData)
}

getData <- function()
{
    library(plyr)
    library(foreign)
    library(tiff)
    library(data.table)
    setwd('/Users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data')
    x <- 1:3
    y <- 0:1
    ImRow <- 0:3
    ImCol <- 0:2
    data <- data.frame()
    countData <- data.frame()
    for(i in x)
    {
        for(j in y)
        {
            countFilename <- paste(getwd(), '/File - Microwell Cell Count Stats 13/x', i, '_y', j, '.arff', sep='')
            print(countFilename)
            countTemp <- read.arff(countFilename)
            countTemp <- reorganizeTable(countTemp)
            countTemp$x <- i
            countTemp$y <- j
            countData <- rbind(countData, countTemp)

            temp <- data.frame()
            for(r in ImRow)
            {
                for(c in ImCol)
                {
                    dataFilename <- paste(getwd(),'/File - Microwell Intensities 13/x', i, '_y', j, '_ImCol', c, '_ImRow', r, '.arff', sep='')
                    print(dataFilename)
                    d <- read.arff(dataFilename)
                    d <- reorganizeTable(d)
                    d$ImRow <- r
                    d$ImCol <- c
                    temp <- rbind(temp, d)
                }
            }
            temp$x <- i
            temp$y <- j
            data <- rbind(data, temp)
        }
    }

    countData$Device <- getDevice(x=countData$x, y=countData$y)
    countData$Image <- getImage(ImRow=countData$ImRow, ImCol=countData$ImCol)
    countData$ID <- getID(device=countData$Device, image=countData$Image, ROI=countData$ROI)
    countData$Single <- countData$Count.Max == 1
    countData$Zero <- countData$Count.Avg == 0
    countData$Virus <- NA
    countData[countData$ID < 37000,]$Virus <- 'N1'
    countData[countData$ID >= 37000,]$Virus <- 'M51R'
    countData <- countData[with(countData, order(ID)),]

    # Read in the image used for illumination correction
    IF <- readTIFF('/users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data/IF.tif')
    meanIF <- mean(IF)
    countData <- data.table(countData, key=c('ID'))
    countData[,IF.Factor:=IF[Y,X]/meanIF,by=ID]

    data$Device <- getDevice(x=data$x, y=data$y)
    data$Image <- getImage(ImRow=data$ImRow, ImCol=data$ImCol)
    data$ID <- getID(device=data$Device, image=data$Image, ROI=data$ROI)
    data$Cell.Count <- NA
    singleIDs <- countData[countData$Single,]$ID
    zeroIDs <- countData[countData$Zero,]$ID
    data[data$ID %in% singleIDs,]$Cell.Count <- 1
    data[data$ID %in% zeroIDs,]$Cell.Count <- 0
    data$Virus <- NA
    data[data$ID < 37000,]$Virus <- 'N1'
    data[data$ID >= 37000,]$Virus <- 'M51R'
    data <- data[with(data, order(ID,Time)),]
    data <- data.table(data, key=c('ID','Time'))

    return(list(data=data, countData=countData))
}

getDevice2 <- function(row, col)
{
    ret <- data.frame()
    for(c in unique(col))
    {
        for(r in unique(row))
        {
            device <- 0
            if(c > 3)
            {
                if(r %% 4 == 0)
                {
                    device <- r %/% 4
                }
                else
                {
                    device <- r %/% 4 + 1
                }
            } else
            {
                if(r %% 4 == 0)
                {
                    device <- r %/% 4 + 3
                }
                else
                {
                    device <- r %/% 4 + 4
                }
            }
            ret <- rbind(ret, data.frame(row=r, col=c, Device=device))
        }
    }
    return(ret)
}

getImage2 <- function(row, col)
{
    ret <- data.frame()
    for(c in unique(col))
    {
        for(r in unique(row))
        {
            newRow <- 0
            if(r %% 4 == 0)
            {
                newRow <- 4
            }
            else
            {
                newRow <- r %% 4
            }

            newCol <- 0
            if(c %% 3 == 0)
            {
                newCol <- 3
            }
            else
            {
                newCol <- c %% 3
            }

            image <- (newCol-1) * 4 + newRow

            ret <- rbind(ret, data.frame(row=r, col=c, ImRow=newRow, ImCol=newCol, Image=image))
        }
    }
    return(ret)
}

getErrors <- function()
{
    errors <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data/ErrorLog2.txt', header=FALSE, sep='\t', col.names=c('row','col','ROI'))
    errors$ROI <- as.numeric(as.character(errors$ROI))
    device <- getDevice2(errors$row, errors$col)
    image <- getImage2(errors$row, errors$col)
    errors <- merge(errors, image, by=c('row','col'))
    errors <- merge(errors, device, by=c('row','col'))
    errors$ID <- getID(errors$Device,errors$Image,errors$ROI)
    errors <- errors[with(errors, order(ID)), ]

    errors$Virus <- NA
    errors[errors$ID < 37000,]$Virus <- 'N1'
    errors[errors$ID >= 37000,]$Virus <- 'M51R'

    return(errors)
}

writeData <- function(data, countData, errors)
{
    library(foreign)
    write.table(data, '/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterDataFile.txt', row.names=FALSE)
    write.table(countData, '/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterCountDataFile.txt', row.names=FALSE)
    write.table(errors, '/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterErrorsFile.txt', row.names=FALSE)
}

readData <- function()
{
    library(foreign)
    print('Getting data')
    data <- data.table(read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterDataFile.txt', header=TRUE))
    print('Getting countData')
    countData <- data.table(read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterCountDataFile.txt', header=TRUE))
    print('Getting errors')
    errors <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/MasterErrorsFile.txt', header=TRUE)
    return(list(data=data, countData=countData, errors=errors))
}

getAndreaData <- function()
{
    library(foreign)
    print('Getting Andrea\'s Data')
    andreaData <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data/AndreaTable_M51R.txt', header=TRUE)
}

rollmin <- function(piece, thresh)
{
    # Apply a rolling min function to the trajectory
    if(length(piece)<4){return(piece)};
    ret<-rollapply(piece,width=4,min,align='left',fill=NA);
    # pad the min-filtered trajectory to be the same length as the original
    ret[is.na(ret)] <- piece[is.na(ret)];
    # Determine the first point above the threshold as the official start of the trajectory
    i <- which(ret > thresh)[1];
    print(piece)
    print(ret)
    print(i)
    if(is.na(i)){return(piece*0);}
    # If we reached the end of the piece (minus the window) then set the whole thing to zero
    if(i <= (length(piece)-3))
    {
        ret[i:length(ret)] <- piece[i:length(ret)];
    } else # Else set the remainder of points to zero
    {
        ret[i:length(ret)] <- 0;
    }
    # Zero all data before that point and return
    if(i > 1){ret[1:i-1] <- 0;}
    return(ret)
}

st <- function(...)
{
    out <- '';
    for(txt in list(...))
    {
        out <- paste(out, as.character(txt), sep='')
    }
    return(out)
}

preprocessVirusData <- function(data, virusType='M51R') # or 'N1'
{
     # Focus on M51R Data
     data <- subset(data, Virus==virusType)

     # Read in the times at which each frame was acquired and fill in the time column with the image frame number (starting at 1)
     timedata <- read.csv('/Users/jaywarrick/Google Drive/SingleCellLatest/Compiled Data/timeData.dat')
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

     #      # Create a quick plot for a sense check
     #      plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$R.Final)))
     #      single[, lines(Time,R.Final, col='red'), by=ID]
     #      zero[, lines(Time,R.Final, col='black'), by=ID]
     #
     #      plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$G.Final)))
     #      single[, lines(Time,G.Final, col='green'), by=ID]
     #      zero[, lines(Time,G.Final, col='black'), by=ID]
     #
     #      plot(c(),c(),xlim=c(0,47),ylim=c(0,max(single$B.Final)))
     #      single[, lines(Time,B.Final, col='blue'), by=ID]
     #      zero[, lines(Time,B.Final, col='black'), by=ID]

     # Determine the threshold based on 0-cell data
     zeroSummary <- zero[, list(R.Mean=mean(R.Final),G.Mean=mean(G.Final),B.Mean=mean(B.Final)), by=list(Device,Image,ID)]
     zeroSummary <- zeroSummary[, list(R.Mean.Median=median(R.Mean), G.Mean.Median=median(G.Mean), B.Mean.Median=median(B.Mean), R.Mean.StdDev=mad(R.Mean), G.Mean.StdDev=mad(G.Mean), B.Mean.StdDev=mad(B.Mean)), by=list(Device,Image)]
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

     #      # Create a quick plot for a sense check
     #      plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$R)), log='y')
     #      single[, lines(Time,R, col=rgb(1,0,0,0.2)), by=ID]
     #
     #      plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$G)), log='y')
     #      single[, lines(Time,G, col=rgb(0,1,0,0.2)), by=ID]
     #
     #      plot(c(),c(),xlim=c(0,47),ylim=c(1,max(single$B)), log='y')
     #      single[, lines(Time,B, col=rgb(0,0,1,0.2)), by=ID]

     return(list(single=single, zero=zero, thresholds=thresholds))
}
