#' Simulate Power for a null hypothesis involving one or more zero-constraints.
#' @description Monte-Carlo-simulate power for a given model to reject a null hypothesis involving one or more zero-constraints.
#' @details 	This is a wrapper function to the generic \code{simpower} function to obtain power for a given model to reject a null hypothesis involving one or more zero-constraints by a Monte Carlo simulation.
#' @param trueModel An OpenMx model representing the population model.
#' @param restrictions A list of parameters to be restricted to zero under H0.
#' @param N Sample size
#' @param repetitions The number of Monte-Carlo-trials.
#' @param keepModels Store the complete intermediate results in form of the fitted models. Only recommended for debugging purposes. 
#'
#' @export

simPowerZeroRestriction<-function(trueModel, restrictions, N, repetitions, keepModels=F, 
                                  steps=NULL, cluster=NULL)
{
  iterate <- steps
  
  if (is.null(iterate)) {
    return(simPowerZeroRestrictionDelegate(trueModel,restrictions,N,
                                           repetitions,keepModels,cluster=cluster))
  } else {
    
    result <- list()
    combn <- expand.grid(iterate)
    for (i in 1:nrow(combn)) {
      for (j in 1:ncol(combn)) {
        trueModel <- omxSetParameters(trueModel,labels=names(combn)[j],values=combn[i,j])
      }
      result[[i]] <- simPowerZeroRestrictionDelegate(trueModel,restrictions,
                                                     N,repetitions,keepModels,cluster=cluster)
    }
    
    ps <- sapply(result, FUN = power)

    combn <- cbind(combn, ps)
        
    return(list(result=result, combn=combn))
  }
  
}
  
simPowerZeroRestrictionDelegate<-function(trueModel, restrictions, N, repetitions, keepModels=F, 
                                    cluster=NULL)
{
  true.model <- trueModel
  keep.models <- keepModels
 # snowfall.cpus <- snowfallCpus
  h1Model <- true.model
  h0Model <- true.model
  
  h0params <- omxGetParameters(model=h0Model)
  filt <- which(  names(h0params) %in% restrictions)
  
  if (any(h0params[filt]==0)) {
    message("Warning. Specified parameters are already zero in the trueModel.")
    #return()
  }
  
  h0Model <- omxSetParameters(model=h0Model,labels=restrictions,free=F,values=0)
  
  if (!is.null(cluster)) {
    cat("Parallel package init!")
#    sfInit(parallel=T, cpus=snowfall.cpus)
    #sfExportAll()
    #sfSource("power.simulationV2.R")
    #sfClusterEval(require("OpenMx"))
    #sfClusterEval(require("sempower"))
    
    snowfall.cpus <- length(cl)
    
    each.rep <- ceiling(repetitions / snowfall.cpus)
    
    # cat("EACH REP", each.rep)
    x <- rep(list(each.rep), snowfall.cpus)
    print(x)

    t1 <- proc.time()
    
    results <- parLapply( X=x, fun=sfWrapSimPower, 
                         h0Model, h1Model, true.model, 
                         N, keep.models)
    t2 <- proc.time()    

    #sfRemoveAll()
    #sfStop()
    
    #
    # merge
    # 
    result <- c()
    result$statistics <- c()
    result$p <- c()
    result$N <- c()
    
    for (i in 1:length(results)) {
      result$p <- c(result$p,      results[[i]]$p)
      result$N <- c(result$N,      results[[i]]$N)
      result$statistics <- c(result$statistics,      results[[i]]$statistics)
    }
    
    result$elapsed <- t2-t1
    
    result$simulation.type <- results[[1]]$simulation.type
    class(result) <- "simPower"
    #sapply( results, function(x) { return(t(x$statistics));})
    
    return(result)
    #return(
    #  result
    #)
    
  } else {
    return(simPower(h0Model, h1Model, true.model, N=N, 
                    repetitions=repetitions,keepModels=keep.models))
    
  }
  
}