getSingleStats <- function(isVerticalThreshold, x, y, cross, thresh)
{
     ret <- hash()
     ret$n <- length(x)
     ret$meanX <- mean(x)
     ret$meanY <- mean(y)
     ret$sdX <- sd(x)
     ret$sdY <- sd(y)

     total <- ret$n

     if(isVerticalThreshold)
     {
          # Define subpopulations
          plus <- (x > (y*cross + thresh))
          minus <- (x <= (y*cross + thresh))

          # Calculate stats for each subpopulation
          ret$'X % +' <- 100 * (sum(plus) / (total));
          ret$'X n +' <- sum(plus)
          ret$'X meanX +' <- mean(x[plus])
          ret$'X sdX +' <- sd(x[plus])
          ret$'X meanY +' <- mean(y[plus])
          ret$'X sdY +' <- sd(y[plus])

          ret$'X % -' <- 100 * (sum(minus) / (total))
          ret$'X n -' <- sum(minus)
          ret$'X meanX -' <- mean(x[minus])
          ret$'X sdX -' <- sd(x[minus])
          ret$'X meanY -' <- mean(y[minus])
          ret$'X sdY -' <- sd(y[minus])
     }
     else
     {
          # Define subpopulations
          plus <- (y > (x*cross + thresh))
          minus <- (y <= (x*cross + thresh))

          # Calculate stats for each subpopulation
          ret$'Y % +' <- 100 * (sum(plus) / (total))
          ret$'Y n +' <- sum(plus)
          ret$'Y meanX +' <- mean(x[plus])
          ret$'Y sdX +' <- sd(x[plus])
          ret$'Y meanY +' <- mean(y[plus])
          ret$'Y sdY +' <- sd(y[plus])

          ret$'Y % -' <- 100 * (sum(minus) / (total))
          ret$'Y n -' <- sum(minus)
          ret$'Y meanX -' <- mean(x[minus])
          ret$'Y sdX -' <- sd(x[minus])
          ret$'Y meanY -' <- mean(y[minus])
          ret$'Y sdY -' <- sd(y[minus])
     }
     return(ret)
}