###################################################
#######  Correlations  ############################
###################################################

vars <- c('NLFit.Alpha.R', 'Time.Delay.R', 'Time.Rise.R', 'Int.Max.R', 'NLFit.Alpha.G', 'Time.Delay.G', 'Time.Rise.G', 'Int.Max.G');
pops <- list(pop_Alpha.Gp, pop_Delay.Gp, pop_MaxInt.Gp, pop_MaxInt.Gp, pop_Alpha.Gp.G, pop_Delay.Gp, pop_MaxInt.Gp, pop_MaxInt.Gp)

Vresults <- data.frame(Parameter=vars, NLFit.Alpha.R=numeric(4), Time.Delay.R=numeric(4), Time.Rise.R=numeric(4), Int.Max.R=numeric(4), NLFit.Alpha.G=numeric(4), Time.Delay.G=numeric(4), Time.Rise.G=numeric(4), Int.Max.G=numeric(4), stringsAsFactors=FALSE)
for(m1 in 1:length(vars))
{
    for(m2 in 1:length(vars))
    {
        m1Name <- vars[m1]
        m2Name <- vars[m2]
        newPop <- merge(pops[m1], pops[m2])
        
        m1IsInt <- length(grep('Int',m1Name,fixed=TRUE)) > 0
        m2IsInt <- length(grep('Int',m2Name,fixed=TRUE)) > 0
        
        d1 <- newPop[,m1Name]
        d2 <- newPop[,m2Name]
        
        if(m1IsInt)
        {
            #             browser()
            d1 <- log10(d1)
        }
        if(m2IsInt)
        {
            #             browser()
            d2 <- log10(d2)
        }
        cor <- cor.test(d1, d2)
        
        Vresults[Vresults$Parameter==m1Name,m2Name] <- cor$estimate
    }
}
data <- Vresults
Vresults
for(m in vars)
{
    Vresults[,m] <- sprintf('%.2f',Vresults[,m])
}
nums <- Vresults
for(r in 1:8)
{
    for(c in 2:9)
    {
        if(abs(as.numeric(Vresults[r,c])) >= 0.25)
        {
            Vresults[r,c] <- st('\\bf ', Vresults[r,c])
        }
        if(r >= c && c > 1)
        {
            Vresults[r,c] <- st('\\gray ',Vresults[r,c])
        }
    }
}
# Vresults <- Vresults[c(1:7),c(1,3:9)]
Vresults
library(xtable)
print(xtable(Vresults), include.rownames=FALSE)

####################################################
#######  Intact vs Lysed Comparison ################
####################################################

wilcox.Results <- data.frame()

pop <- pop_Alpha.All
p <- 'NLFit.Alpha.R'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_Alpha.Gp.G
p <- 'NLFit.Alpha.G'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_Delay.All
p <- 'Time.Delay.R'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_Delay.Gp
p <- 'Time.Delay.G'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_MaxInt.All
p <- 'Time.Rise.R'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_MaxInt.Gp
p <- 'Time.Rise.G'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_MaxInt.All
p <- 'Int.Max.R'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))

pop <- pop_MaxInt.Gp
p <- 'Int.Max.G'
x1 <- pop[pop$Death == 1,p]
x2 <- pop[pop$Death == 2,p]
mean(x1, ra.rm=TRUE)
mean(x2, ra.rm=TRUE)
t <- wilcox.test(x1,x2, alternative='two.sided')
wilcox.Results <- rbind(wilcox.Results, data.frame(parameter=p, intact=mean(x1), lysed=mean(x2), p.value=t$p.value, alternative=t$alternative))


library(xtable)
fileConn<-file('/Users/jaywarrick/GoogleDrive/SingleCell/LysisTable.tex')
writeLines(print(xtable(wilcox.Results, display=rep('fg', ncol(wilcox.Results) + 1), digits=2), include.rownames=FALSE), fileConn)
close(fileConn)

print(xtable(wilcox.Results, display=rep('fg', ncol(wilcox.Results) + 1), digits=2), include.rownames=FALSE)


##########################################
#######  GFP+ vs GFP- Comparison #########
##########################################

# Max.R vs Alpha.R (Intersect MaxInt and Alpha pops)
delayP <- pop_Delay.Gp$Time.Delay.R
riseP <- pop_MaxInt.Gp$Time.Rise.R
alphaP <- pop_Alpha.Gp$NLFit.Alpha.R
intP <- pop_MaxInt.Gp$Int.Max.R
pop <- pop_MaxInt.Gp
lysedP <- length(which(pop$Death==2)) # Use MaxInt population because contains only cells where either a cell plateaued or lysed
totP <- nrow(pop)

delayN <- pop_Delay.Gn$Time.Delay.R
riseN <- pop_MaxInt.Gn$Time.Rise.R
alphaN <- pop_Alpha.Gn$NLFit.Alpha.R
intN <- pop_MaxInt.Gn$Int.Max.R
pop <- pop_MaxInt.Gn
lysedN <- length(which(pop$Death==2)) # Use MaxInt population because contains only cells where either a cell plateaued or lysed
totN <- nrow(pop)

p.delay <- wilcox.test(delayP,delayN)
p.rise <- wilcox.test(riseP,riseN)
p.alpha <- wilcox.test(alphaP,alphaN)
p.int <- wilcox.test(intP,intN)
p.lysed <- prop.test(c(lysedP,lysedN),c(totP,totN))

comparisonSummary <- data.frame(Parameter=c('Delay-Time','Rise-Time','Prod.-Rate','Max-Int.','% Lysed'), 'GFP.Pos'=c(mean(delayP),mean(riseP),mean(alphaP),mean(intP),(100*lysedP/totP)), 'GFP.Neg'=c(mean(delayN),mean(riseN),mean(alphaN),mean(intN),(100*lysedN/totN)), 'p.value'=(c(p.delay$p.value,p.rise$p.value,p.alpha$p.value,p.int$p.value,p.lysed$p.value)))

library(xtable)
fileConn<-file('/Users/jaywarrick/GoogleDrive/SingleCell/GFPTable.tex')
writeLines(print(xtable(comparisonSummary, display=c(rep('fg', ncol(comparisonSummary)),'e'), digits=1, include.rownames=FALSE), fileConn))
close(fileConn)
print(xtable(comparisonSummary, display=c(rep('fg', ncol(comparisonSummary)),'e'), digits=3), include.rownames=FALSE)

rownames(data) <- data$Parameter
data2 <- data[, 2:length(data[1,])]
data_mat <- as.matrix(data2, rownames.force=TRUE)

rownames(nums) <- nums$Parameter
nums2 <- nums[, 2:length(nums[1,])]
nums_mat <- as.matrix(nums2, rownames.force=TRUE)

grays <- function (n, alpha = 1) 
{
    alpha <- 1
    if ((n <- as.integer(n[1L])) > 0L) {
        even.n <- n%%2L == 0L
        k <- n%/%2L
        l1 <- k + 1L - even.n
        l2 <- n - k + even.n
        lower <- round(seq.int(0.5, ifelse(even.n, 1/k, 0), length.out = l1),2)
        upper <- round(seq.int(0, 0.5, length.out = l2)[-1L],2)
        ret = c(if (l1 > 0L) rgb(r = 1-lower, g = 1-lower, b = 1, alpha = alpha), 
                if (l2 > 1) rgb(r = 1, g = 1-upper+min(upper), b = 1-upper+min(upper), alpha = alpha))
        ret[length(ret)] <- gray(0);
        print(lower)
        print(upper)
        return(ret)
    }
    else character()
}

getColors <- function(n, alpha = 1, unityColor=rgb(0.5,0.5,0.5,1)) 
{
    if(n%%2 == 0)
    {
        n <- n + 1
    }
    halfn <- (n-1)/2
    
    upper <- (0:halfn)/halfn
    lower <- rev(upper)
    print(lower)
    print(upper)
    ret <- c(rgb(r=1-lower, g=1-lower, b=1), rgb(r=1, g=1-upper, b=1-upper)[-1])
    ret[length(ret)] <- unityColor
    return(ret)
}

getBreaks <- function(n)
{
    if(n%%2 == 0)
    {
        n <- n + 1
    }
    ret <- 2*(0:(n))/(n) - 1
    return(ret)
}

library(pheatmap)
pdf('/Users/jaywarrick/GoogleDrive/SingleCell/Figures/Heatmap.pdf', width=4.5, height=4)
temp <- data_mat[,c(8,7,6,5,4,3,2,1)]
for(i in 1:8)
{
    for(j in 1:8)
    {
        if(i + j > 9)
        {
            temp[i,j] <- NA
        }
    }
}
pheatmap(temp, number_format = "%.2f", show_rownames=F, show_colnames=F, border_color='white', fontsize=10, fontsize_number=8, cluster_rows=F, cluster_cols=F, color=getColors(100), breaks=getBreaks(100), display_numbers=T)
dev.off()

