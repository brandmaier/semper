
modelPredictedCov <- function(model, type="raw") {
  run <- model
  fit <- NULL
  data <- NULL
  manifests <- run@manifestVars

  
  if (type=="raw") {
    fake.data <- data.frame(matrix(1:(length(manifests) * 2), 
                                   ncol = length(manifests)))
    names(fake.data) <- manifests
    run@data <- mxData(fake.data, type = "raw")
  } else if (type=="cov") {
    fake.data <- diag(nrow = length(manifests))
    colnames(fake.data) <- manifests
    rownames(fake.data) <- manifests
    run@data <- mxData(fake.data,type="cov",numObs = 1000)  
  } else {
    stop("Unknwon type.")
    return(NA);
  }
  
  fit <- mxRun(run, useOptimizer = F, silent = T)
  cov <- attr(fit$output$algebras[[1]], "expCov")  
  
  return(cov);
}

