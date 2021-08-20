
# matrix trace
mtr <- function(m) {sum(diag(m))}

#
# helper for numeric computation
#
loss.func <- function(x, chiByN, N) {
  ratio <- 1/(1-x)
  loss = (chiByN-(ratio-log(ratio)-1))^2
  return(loss)
}

compute.numeric <- function(chiByN) {
  ecr <- optimise(loss.func, c(0,1), chiByN=chiByN, tol = 1e-8)
  return(ecr$minimum)
}

#
#' generic Effective Curve Reliability (ECR)
#'
#' @param model A semper model
#' @param constraints A list of parameter names assumed to be zero in the null hypothesis
#'
#' @export

ecr.generic <- function(model, constraints="slope.variance") {
  
  dataMat <- getExpectedCovariance(model)
  
  model0 <- model

  if (inherits(model,"semper")) {  
  
    if (!all(constraints %in% attributes(model0)$names)) stop("Not all constraints in model!")
    
  if (is.vector(constraints)) {
    for (cn in constraints)
      model0[cn] <- 0
  } else {
    model0[constraints] <- 0
  }
    
  }  else if (inherits(model,"MxRAMModel")) {
    model0 <- omxSetParameters(model, labels=constraints, values=0)
  } else {
    stop("Unknown model type")
  }
  
  dataMat0 <- getExpectedCovariance(model0)
  
  d0 <- dataMat0
  d1 <- dataMat
  
  if (all(d1==d0)) stop("Error! Expected covariance matrices are identical!")
 
  
  
  # do not multiply by N because this is chi^2/N
  dchi <- (log(det(d0))+mtr(solve(d0)%*%d1)-log(det(d1))-dim(d0)[1])
  
  ecr.new <- compute.numeric(dchi)
  
  return(ecr.new)
  
}

#
# turn ECR back into an effective erro
# assuming that effect size is in absolute terms!
# 
effective.error.from.ecr<-function(abs.effect, ecr.value)
{
  abs.effect <- abs(abs.effect)
  return( abs.effect*(1/ecr.value-1))
}

#'
#' @export

effectsize <- function(model, constraints="slope.variance") {
  
  if(is.null(model)) stop("No model given!")
  
  dataMat <- getExpectedCovariance(model)
  
  model0 <- model
  if (!all(constraints %in% attributes(model0)$names)) stop("Not all constraints in model!")
  
  if (is.vector(constraints)) {
    for (cn in constraints)
      model0[cn] <- 0
  } else {
    model0[constraints] <- 0
  }
  
  dataMat0 <- getExpectedCovariance(model0)
  
  d0 <- dataMat0
  d1 <- dataMat
  
  if (all(d1==d0)) stop("Error! Expected covariance matrices are identical!")
  
  #dchi <- log(det(d0))+mtr(solve(d0)%*%d1)+log(2*pi)
  #dchi <- log(det(d1))+mtr(solve(d1)%*%d0)
  
  
  # do not multiply by N because this is chi^2/N
  dchi <- (log(det(d0))+mtr(solve(d0)%*%d1)-log(det(d1))-dim(d0)[1])
  
#  dn <- 1/(1+1/dchi)
  
  #dn <- - dchi / (dchi - 1)
  
  return(dchi)
}
