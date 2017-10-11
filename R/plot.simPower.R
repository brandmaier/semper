
plot.simPower <- function(psim, prng=seq(0.4,0.95,0.05), target.power=NULL, lw=2, lty=1,
                           xlab="Sample Size", ylab="Statistical Power",
                           main="Monte Carlo Power Simulation", add=F, ...)
{
  
  if (psim$simulation.type=="power") {
    
    pow <- power(psim)
    hh <- hist(psim$p)
    brk <- c(hh$breaks,0.05)
    hist(psim$p, breaks=brk,main=paste("Monte Carlo p-values (Estimated Power=",pow,")",sep=""),xlab="p")
    abline(v=0.05,lty=2)
    return();
  }
 
  N.hat <- sampleSize(psim, power=prng)
  if (!add)
  plot(N.hat,prng,type="n", xlab=xlab,
       ylab=ylab,main=main,
       frame.plot=F, ...)
  lines(N.hat,prng, lw=lw,lty=lty)
  
  if (!is.null(target.power)) {
   target.N <- N.hat[prng==target.power]
   abline(v=target.N,lty=2)
   abline(h=target.power, lty=2)
   text( x=N.hat[prng==target.power], y = prng[1],labels=paste("N>=",ceiling(target.N)))
  }
} 