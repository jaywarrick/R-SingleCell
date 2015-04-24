getDensityColors <- function(x, y)
{
#     library(hexbin)
    xunique <- unique(x)
    yunique <- unique(y)
    xync <- data.frame(x=c(), y=c(), n=c(), c=c())
    for(xu in xunique)
    {
        for(yu in yunique)
        {
            i <- sum(x==xu & y==yu)
            xync <- rbind(xync, data.frame(x=xu, y=yu, n=i))
        }
    }
    xync <- subset(xync, n > 0)
    theColors <- rev(rainbow(max(xync$n), end=0.6))
#     theColors <- plinrain(5, beg=50, end=200)
#     plot(c(), c(), xlim=c(0, length(theColors) + 1), ylim=c(0, length(theColors) + 1))
    for(i in 1:length(theColors))
    {
        points(i, i, pch=21, col=rgb(0,0,0,0), bg=theColors[i])
    }
    xync$c <- theColors[xync$n]
    return(list(colors=theColors, data=xync))
}