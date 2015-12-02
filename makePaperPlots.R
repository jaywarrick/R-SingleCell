rm(list=ls())
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/plotHelperFunctions.R')
library(Hmisc)

# Read the raw and summary data
data <- read.table('/users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Data_1Cell.txt', header=TRUE)
params <- read.table('/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/TrajectoryParameters_20131027.txt', header=TRUE)
thresholds <- read.table('/Users/jaywarrick/Google Drive/SingleCellLatest/Processed Data/M51R_Thresholds.txt', header=TRUE)

# Define all the populations that will be potted
pop_MaxInt.Gp <- subset(params, Flag.R & Flag.G & (Death==1 | Death==2)) # Remove no plateau rois
pop_MaxInt.Gn <- subset(params, Flag.R & !Flag.G & (Death==1 | Death==2)) # Remove no plateau rois
pop_MaxInt.All <- rbind(pop_MaxInt.Gp, pop_MaxInt.Gn)

pop_Delay.Gp <- subset(params, Flag.R & Flag.G)
pop_Delay.Gn <- subset(params, Flag.R & !Flag.G)
pop_Delay.All <- rbind(pop_Delay.Gp, pop_Delay.Gn)

# For Red Fits Only
pop_Alpha.Gp <- subset(params, Flag.R & Flag.G & NLFit.Flag.Quality.R)
pop_Alpha.Gn <- subset(params, Flag.R & !Flag.G & NLFit.Flag.Quality.R)
pop_Alpha.All <- rbind(pop_Alpha.Gp, pop_Alpha.Gn)

# For Green Fits Only
pop_Alpha.Gp.G <- subset(params, Flag.R & Flag.G & NLFit.Flag.Quality.G) # This is the only population of the three that really makes sense to use
pop_Alpha.Gn.G <- subset(params, Flag.R & !Flag.G & NLFit.Flag.Quality.G)
pop_Alpha.All.G <- rbind(pop_Alpha.Gp.G, pop_Alpha.Gn.G)

# For Red & Green Fits Only
pop_Alpha.Gp.RG <- subset(params, Flag.R & Flag.G & NLFit.Flag.Quality.R & NLFit.Flag.Quality.G)  # This is the only population of the three that really makes sense to use
pop_Alpha.Gn.RG <- subset(params, Flag.R & !Flag.G & NLFit.Flag.Quality.R & NLFit.Flag.Quality.G)
pop_Alpha.All.RG <- rbind(pop_Alpha.Gp.RG, pop_Alpha.Gn.RG)

# In general, look at how good the fits were
sum(params[params$Flag.R,]$NLFit.R2.R>=0.9, na.rm=TRUE)/length(unique(data[data$Flag.R,]$ID))
sum(params[params$Flag.R & params$Flag.G,]$NLFit.R2.G>=0.9, na.rm=TRUE)/length(unique(data[data$Flag.R & data$Flag.G,]$ID))
mean(params[params$Flag.R,]$NLFit.R2.R, na.rm=TRUE)
mean(params[params$Flag.R & params$Flag.G,]$NLFit.R2.G, na.rm=TRUE)
sum(params$NLFit.R2.R>=0.9, na.rm=TRUE)/length(unique(data$ID))
mean(params$NLFit.R2.R, na.rm=TRUE)
mean(params$NLFit.R2.G, na.rm=TRUE)

# plotting parameters
G <- 'darkgreen'
R <- 'darkred'
H <- 4.7
W <- 4.5

source('/Users/jaywarrick/Public/DropBox/R-SingleCell/plotScriptFACS.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/plotScriptMicroscope.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/plotScriptMicrowell.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_6.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_7.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_8.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_9.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_10.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_11.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Figure_12.R')
source('/Users/jaywarrick/Public/DropBox/R-SingleCell/Tables.R')

