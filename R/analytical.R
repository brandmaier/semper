#effectiveError <- function(ti, err, icept, cov)
#{
#	
#M <- length(ti)
#eff.err <- (cov**2*(sum(ti)**2-M*sum(ti*ti))+(2*cov*sum(ti)+M*icept)*err+err**2) / (sum(ti*ti)*(M*icept+err)-icept*sum(ti)**2)#
#
#return(eff.err)
#}


powerFromError <- function(lgcm, alpha=.05){

slope <- lgcm$slope.variance
eff.err <- effective.error(lgcm)
N <- lgcm$sample.size
  
crit <- qchisq(1-alpha,df=1,ncp=0)	
lambda <- log(slope+eff.err)-log(eff.err)+eff.err/(slope+eff.err)-1
lambda <- (lambda*N)

# Area under non-central chi^2 
power = 1-pchisq(crit, df=1, ncp=lambda)

return(power)
}

powerFromError2df <- function(lgcm, alpha=.05) {
  
  ncp <- ncp2df(lgcm)
  
  # Area under non-central chi^2 
  crit <- qchisq(1-alpha,df=2,ncp=0)	# 5.99 for alpha=.05
  power = 1-pchisq(crit, df=2, ncp=lambda)
  
  return(power)
}

ncp2df <- function(lgcm)
{
  ti<- lgcm$timepoints
  err <- lgcm$residual.variance
  icept <- lgcm$intercept.variance
  slope <- lgcm$slope.variance
  rho <- lgcm$intercept.slope.covariance
  N <- lgcm$sample.size 
  
  tisum <- sum(ti)
  tisqsum <- sum(ti*ti)
  numT <- length(ti)

  lambda <- tisum/tisqsum
  e1 = err / (numT-tisum*tisum/tisqsum)
  e2 = err/tisqsum
  
  Sigma = matrix(
    c(
      icept+e1, lambda*icept+rho,
      lambda*icept+rho, lambda^2*icept+2*lambda*rho+e2+slope
      )
    ,nrow=2,byrow=T
    )
  
  SigmaRes = matrix(
    c(
      icept+e1, lambda*icept,
      lambda*icept, lambda^2*icept+e2
    )
    ,nrow=2,byrow=T
  ) 
  
  ncp <-  -log(det(Sigma)) + log(det(SigmaRes)) + 
                   sum(diag( solve(SigmaRes) %*% Sigma )) -2 
  
  ncp <- ncp*N
  
  return(ncp)
}

