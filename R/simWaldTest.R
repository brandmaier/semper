#' Simulate power of the Wald test
#' @param lmodel a \code{semper} model
#' @param repetitions Number of Monte Carlo trials
#' @param A parameter name assumed to be zero in the null hypothesis
#' @param sample.size The sample size 
#' @export

simPowerWaldTest <- function(lmodel, repetitions=1000, paramName="slopevariance", 
                             sample.size=NULL, simfunc=simulateData) {
  
  cnt <- 0
  model <- getOpenMxRepresentation(lmodel)
  true.model <- model
  zs <- rep(NA, repetitions)
  ests <- rep(NA, repetitions)
  Ns <- rep(NA, repetitions)
  simulation.type <- "power"
  ps <- rep(NA, repetitions)
  for (i in 1:repetitions) {

   if (is.null(sample.size)) {
    N <- lmodel$sample.size
   } else {
	if (length(sample.size)==1) {
		N <- sample.size
	} else if (length(sample.size)==2) {
		N <- floor(runif(1, sample.size[1],sample.size[2]))
		simulation.type <- "samplesize"
	} else {
		stop("Choice of N is Not supported!")
	}
   }

   simdata <- simfunc(true.model,N )

   model <- mxModel(true.model, mxData(simdata,type="raw"))
   run <- mxRun(model,silent = TRUE) 
   sm <- summary(run)
   lzs <- sm$parameters$Estimate / sm$parameters$Std.Error
   pid <- sm$parameters$name==paramName
   z <- lzs[pid]
   if (abs(z) > 1.96) cnt<-cnt+1
    zs[i] <- z
   ests[i] <- sm$parameters$Estimate[pid] 
   Ns[i] <- N
   ps[i] <- pnorm(q = z,lower.tail = FALSE)
  }
  xx<- list(z=zs,estimates=ests,N=Ns,simulation.type=simulation.type,p=ps, statistic.type="z")
  class(xx) <- "simPower"
	return(xx)
}
