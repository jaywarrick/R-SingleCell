logicle <- function(x, transitionPoint, linearUnitsPerOrder)
{
    linearDifference = x - transitionPoint;
    valsToAdjust <- linearDifference <= 0;
    ordersDifferenceOnDisplay = linearDifference / linearUnitsPerOrder;
    linearDifference[valsToAdjust] <- transitionPoint * 10^(ordersDifferenceOnDisplay[valsToAdjust]);
    x[valsToAdjust] <- transitionPoint*10^(-1*((transitionPoint-x[valsToAdjust])/linearUnitsPerOrder))
    
    return(x)
}