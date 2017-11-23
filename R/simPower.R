#'
#' Monte-Carlo power simulation
#'
#' Power simulation of a likelihood ratio test between two
#' models given simulated data from a true model in R/OpenMx
#' Assume you have a H1 model that test whether a pretest and a posttest score
# mean are different in a sample. The corresponding H0 model states that
# there is no difference between the scores, and therefore constitutes a
# restricted and nested version of the H1 model. Given a true model that
# states the difference in the scores, data is repeatedely simulated
# and the power to detect this difference can be calculated given a
# chosen sample size. Alternatively, the required sample size can be
# calculated given a chosen power.
#'
#'
#' @param h0Model Model representing the null hypothesis
#' @param h1Model Model representing the alternative hypothesis
#' @param populationModel}{True population model used to generate data from.
#' @param N Sample size. For RAM models, this is either (1) a scalar representing a fixed sample size, (2) a list with two items representing a minimum and a maximum for randomly drawn sample sizes (when sample size is the dependent variable), or (3) for multiple group models, this is a named list of sample sizes for the sub models where the names must match the submodels' names.
#'
#'@parampsim SimPower result object
#'@paramprng Range for plotting p values
#'@paramtarget.power If not null, draw a line to indicate sample size for a specified target power in a power curve plot.
#'@paramlw line width for power curve
#'@paramlty line type for power curve
#'
#' @seealso \code{\link{simPowerZeroRestriction}}
#'
#' @aliases simPower print.simPower plot.simPower
#'
#' @examples
#' simPower(h0Model,h1Model, populationModel, N=c(20,100),repetitions=100, keepModels=F)
#'
#' plot(psim, prng=seq(0.4,0.95,0.05), target.power=NULL, lw=2, lty=1,
#'                          xlab="sample size", ylab="statistical power",
#'                           main="Monte Carlo Power Simulation", add=F,...)
#'						   
#' print(psim, ...)	
#'
#' @author: Andreas Brandmaier
#'
#' @detail:  Given a true model, data is repeatedely simulated
#' and the power to detect this difference is calculated for a
#' chosen sample size. Alternatively, the required sample size can be
#' calculated for a chosen power.
#' A simplified wrapper for the simulation of zero constraints is \code{\link{simPowerZeroRestriction}}.
#'
#' @export



simPower <- function(h0Model,h1Model, 
	populationModel, N=c(20,100),repetitions=100, keepModels=FALSE,
	simfunc=simulateData, simargs=NULL)
{
  #stopifnot(is.integer(repetitions))
  N.range <- N

  
	simulation.type = "default"
  
	with.submodels <-  length(populationModel@submodels)>0;
	
	N.names <- NULL
	
  #cat("Model has submodels ",with.submodels,"\n")
	
	if (with.submodels) {
		
		if (!is.list(N.range)) {
			warning("N.range must be a named list!")
			return(null);
		}
		
		N.names <- names(N.range)
	
		# sanity check: do the sample size names match the population sub models
		if (0!=length(setdiff(N.names, names(populationModel@submodels)))) {
			warning("ERROR! N.range must be a list of names matching the submodels!");
			return(null);
		}
		  N.samples <- list()
		  tN.range <- N.range[[1]]
		  
		  
		  
		  if (length(tN.range)==1) {
		    simulation.type="power"
		    N.samples <- list(N.range)
		    for (i in 2:repetitions) {
		      N.samples <- append(N.samples, list(N.range))
		    }
		    
		  } else  if (length(tN.range)==2) {
		    simulation.type="samplesize"
		    N.samples <- list()
		    
		  #  if (repetitions>1)
		    for (i in 1:repetitions) {
		      N.rand <- list()
		      for (name in N.names)
		      {
		        N.rand[[name]] <- floor(runif(1, min=N.range[[name]][1], max=N.range[[name]][2]+1))
		      }
		      N.samples <- append(N.samples, list(N.rand))
		    }
		    
		  } else {
		    stop("Unknown or unsupported parameterization of N");
		  }
		
		
	} else {    # -- without submodels ------
	 if (length(N.range)==2) {
	   simulation.type="samplesize"
		N.samples <- floor(runif(repetitions, min=N.range[1], max=N.range[2]+1))
	 } else if (length(N.range==1)) {
	   simulation.type="power"
		N.samples <- rep(N.range, repetitions)
	 }
	}
	
	vals <- rep(NA, repetitions)
	p <- rep(NA, repetitions)
	# simulate data
	#for (i in 1:repetitions) {
	#  innerPower()
	#}
  t1 <- proc.time()
  result <- sapply(FUN=innerPower,X=N.samples, populationModel,
                   h0Model=h0Model,h1Model=h1Model, N.names=N.names,
                   with.submodels=with.submodels, keep.models=keepModels,
				   simfunc=simfunc, simargs=simargs,
                   simplify=F)
  t2 <- proc.time()
  
  #if (keep.models) {
#  	models <- result$models
#  	result <- result$result	
#  }
# return(result)

  p <- sapply(result, FUN=function(x){return(x$p)}, simplify=T)
  vals <- sapply(result, FUN=function(x){return(x$diffLL)}, simplify=T)
	
#  return(result)
  ret <- list(statistics=vals, N=N.samples, p =p, 
	  simulation.type=simulation.type, statistic.type="chisq")
  
  if (keepModels)
  {
    ret$models1 <- sapply(result, FUN=function(x){return(x$model1)}, simplify=T)
    ret$models2 <- sapply(result, FUN=function(x){return(x$model2)}, simplify=T)
    
      	#ret$models <- result$models
  }

  ret$elapsed <- t2-t1
  
  class(ret) <- "simPower"
  return(ret)
}

innerPower <- function(N.samples, populationModel, h0Model,h1Model, N.names, with.submodels, keep.models=F,
	simfunc=simulateData, simargs=NULL)
{
	try({
  if (with.submodels) {
    model1 <- mxModel(h0Model)
    model2 <- mxModel(h1Model)
    
    for (name in N.names) {
      #browser()
	  if (!is.null(simargs)) {
      	data <- simfunc(populationModel[[name]], N.samples[[name]])
  	  } else {
        data <- simfunc(populationModel[[name]], N.samples[[name]], simargs)  	  	
  	  } 
      #browser()
      model1[[name]] <- mxModel(h0Model[[name]], data=mxData(data,type="raw"))
      model2[[name]] <- mxModel(h1Model[[name]], data=mxData(data,type="raw"))
      
    }
  } else {
    data <- simulateData(populationModel, N.samples)
    model1 <- mxModel(h0Model, data=mxData(data,type="raw"))
    model2 <- mxModel(h1Model, data=mxData(data,type="raw"))
  }
  
  
  run1 <- mxRun(model1,silent=T)
  run2 <- mxRun(model2,silent=T)
  
  
  #print(mxCompare(run1,run2))
  #vals[i] <-  mxCompare(run1,run2)$diffLL[2]
  #p[i] <- mxCompare(run1,run2)$p[2]
  
  #if (i%%10==0) {cat("Finished",i," runs\n")}
 # llr = run1$objective@result-run2$objective@result # OMX 1.0
  llr = run1$output$Minus2LogLikelihood-run2$output$Minus2LogLikelihood
  df = length(omxGetParameters(run2))-length(omxGetParameters(run1))
  p = 1-pchisq(llr,df)	
  
  if (keep.models) {
    return(list(p=p,diffLL=llr,model1=run1, model2=run2));
  } else {
  return(list(p=p,diffLL=llr));
  }
  
  })
  
  return(list(p=NA,diffLL=NA));
 

}




