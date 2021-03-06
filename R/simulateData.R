#' Simulate Data
#' @description Simulating data from a SEM by drawing from the multivariate normal as described by the model covariance matrix and means vector.
#'
#' 
#' @examples
#' simulateData(
#'	toOpenMx(lgcm(slope.variance=1, 
#'	intercept.variance=2, timepoints=0:8), N=100))
#'
#'@param model An OpenMx model.
#'@param N Desired sample size.
#'@param format Either "long" or "wide".
#' @returns A \code{data.frame} with columns relating to the observed variables of the \code{model} and rows to simulated observations.
#' @author Andreas M. Brandmaier
#'
#' @export

simulateData <- function(model, N, format="wide", seed=NULL) {
  
if (inherits(model, "semper")) {
  model <- toOpenMx(model)
}
  
run <- model
fit <- NULL
data <- NULL
manifests <- run@manifestVars

fake.data <- data.frame(matrix(1:(length(manifests)*2),ncol=length(manifests)))
names(fake.data)<- manifests
run@data <- mxData(fake.data,type="raw")

fit <- mxRun(run,useOptimizer=F,silent=T)
#cov <- fit$objective@info$expCov    #omx 1.0
#mean <- fit$objective@info$expMean  #omx 1.0
cov <- attr(fit$output$algebras[[1]],"expCov")
mean <- attr(fit$output$algebras[[1]],"expMean")

if (!is.null(seed)) {
	data <- R.utils::withSeed(mvrnorm(n=N, mu=mean, Sigma=cov),seed=seed)
} else {
	data <- mvrnorm(n=N, mu=mean, Sigma=cov)
}

dimnames(data)[2] <-  dimnames(fit@data@observed)[2]

if (format=="wide") {
  # it's all OK
} else if (format=="long") {
  M <- length(manifests)
  lngdat <- as.vector(t(data))
  data <- data.frame(id=as.factor(rep(1:N, each=M)), obs=as.factor(rep(1:M,N)), x=lngdat)
} else {
  stop("Unknown format")
}

return(data)
}