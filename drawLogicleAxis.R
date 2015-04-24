drawLogicleAxis <- function(ticks=c(-10,1,10,100,1000,10000,100000,1000000), tickLabels=ticks, axisNum=1, transition=1, linLogRatio=1)
{
    #tickLabels <- tickLabels[1:length(ticks)]
    # Create the tick labels and tick locations
#     for(i in 1:length(ticks))
#     {
#         if(ticks[i] == transition)
#         {
#             tickLabels[i] = paste(as.character(tickLabels[i]), '*', sep ='')
#         }
#     }
    ticks <- logicle(ticks, transition, linLogRatio)
    if(axisNum == 2)
    {
        axis(axisNum, at=ticks, labels=tickLabels, las=2)
    }
    else
    {
        axis(axisNum, at=ticks, labels=tickLabels)
    }
    
}