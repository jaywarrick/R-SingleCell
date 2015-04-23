minor.ticks.axis <- function(ax,n,t.ratio=0.5,at=c(),labels=c(),...){

     lims <- par("usr")
     if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]

     if(length(at) == 0)
     {
          major.ticks <- pretty(lims,n=5)
          if(missing(mn)) mn <- min(major.ticks)
          if(missing(mx)) mx <- max(major.ticks)

          major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
     }
     else
     {
          major.ticks <- at
     }

     if(length(labels)==0)
     {
          labels <- sapply(major.ticks,function(i)
               as.expression(bquote(10^ .(i)))
          )
     }

     axis(ax,at=major.ticks,labels=labels,...)

     n <- n+2
     minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
     minors <- minors[-c(1,n)] #grab everything but the first and last element

     minor.ticks = c(outer(minors,major.ticks,`+`))
     #minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


     axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

n_sig <- round(pi*8^2)
n_rad <- pi*11^2
n_rect <- 33*33
n_tot <- seq(n_sig+1, n_rect)

alpha <- sqrt((18000/0.6)/16384) # estimate for alpha for our camera
beta <- 20 # sqrt(n_rad)*1 # estimate of pixel-to-pixel bg noise (i.e., beta) determined by observing that measurements for n_rad pixels of pure background exhibit noise on the order of 1 unit.
I_sig <- 100000
I_bg <- 0
n_bg <- n_tot-n_sig

SNR <- ((n_sig*I_sig - n_tot*I_bg)/(n_tot))/(sqrt((alpha*sqrt(I_sig)/sqrt(n_sig) + beta/sqrt(n_bg))^2 + ((beta^2)*n_bg)))
plot(n_sig/n_tot, log10(SNR/max(SNR)), yaxt='n', ylab="SNR/max(SNR)", xlab=expression('n'['sig']*'/n'['tot']), type='l', lwd=2)
minor.ticks.axis(2,9, at=log10(c(.001,.01,.1,1)), labels=c(.001,.01,.1,1))

#Plot the point corresponding to a single cell in a microwell
i1 <- which.min((n_sig/n_tot - n_sig/(n_rect))^2)
points(n_sig/n_tot[i1],log10(SNR/max(SNR))[i1], pch=21, col='black', bg='black')
# Plot the point corresponding to a single cell in a radial ROI
i2 <- which.min((n_sig/n_tot - n_sig/(n_rad))^2)
points(n_sig/n_tot[i2],log10(SNR/max(SNR))[i2], pch=21, col='black', bg='white')

SNR[i2]/SNR[i1]
