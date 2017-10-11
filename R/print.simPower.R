
print.simPower <- function(x, ...)
{
	
	
  if (x$simulation.type=="power") {
    pow <- sum(x$p<0.05)/length(x$p)
    err <- sqrt(pow*(1-pow)/length(x$p))
   cat("Estimated power at alpha=5%: ",pow," (Std. error ",err,")\n"); 
  } else {
   cat("Estimated samplesize for power=.8 ", 
   ceiling(sampleSize(x, power=.8)),"\n"); 
  }
  cat("Based on ",length(x$p)," Monte Carlo draws\n");
}
