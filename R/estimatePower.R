

power <- function(psim, alpha=.05, exclude=c())
{
  selector <- rep(TRUE,length(psim$p))
  if ("negstat" %in% exclude) {
	  selector <- psim$statistics>=0
  }

  success <- (psim$p[selector] <alpha)
  
  return( sum(success)/sum(selector))
}
