#'
#' @export

sampleSize <- function(psim, power=0.8, alpha=0.05)
{
  
  stopifnot(class(psim)=="simPower")
  
  stopifnot(length(unique(psim$N))>1)
  
  temp <- list()

  temp$success <- (psim$p < 0.05)

  if (psim$statistic.type=="chisq") {
  	temp$success <- (psim$p < alpha)
  } else if (psim$statistic.type=="z") {
	crit <- qnorm(1-alpha)
	temp$success <- (psim$statistics > crit)
  } else {
   stop("Unknown statistic")
  }
  
  temp$N <- psim$N
  myglm <- glm( success~N, data=temp, family=binomial("logit"))
  smyglm <- summary(myglm)
  
  N.estimate <- (log(power/(1-power))-smyglm$coefficients[1,1])/smyglm$coefficients[2,1]
  
  #power <- seq(0.1,0.9,0.01)
  #N.hat <- (log(power/(1-power))-smyglm$coefficients[1,1])/smyglm$coefficients[2,1]
  
  return(N.estimate)
  
}