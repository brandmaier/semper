#' Create a latent growth curve model
#' @description This function creates a lgcm object representing a linear latent growth curve model
#' @param timepoints Vector of time points of measurement occasions
#' @param sample.size Number of cases
#' @param residual.variance Residual variance or inverse of the precision of measurement
#' @param slope.variance Individual differences in slopes
#' @param intercept.variance Individual differences in intercept

lgcm <- function(timepoints=0, intercept.variance=0, slope.variance=0,
                 residual.variance=0, 
                 intercept.slope.covariance = 0, sample.size = 0,
                 intercept.mean = 0, slope.mean = 0)
{
  # set members
  lgcm <- list()
  lgcm$timepoints <- sort(timepoints)
  lgcm$intercept.variance <- intercept.variance
  lgcm$slope.variance <- slope.variance
  lgcm$residual.variance <- residual.variance
  lgcm$total.study.time <- lgcm$timepoints[length(lgcm$timepoints)]
  lgcm$intercept.slope.covariance <- intercept.slope.covariance
  lgcm$sample.size <- sample.size
  
  lgcm$intercept.mean <- intercept.mean
  lgcm$slope.mean <- slope.mean
  
  # some calculations
  lgcm$sumti <- sum(lgcm$timepoints)
  lgcm$sumtisq <- sum(lgcm$timepoints^2)
  lgcm$num.timepoints <- length(timepoints)
  
  # set return class
  class(lgcm) <- c("semper","lgcm")

  # return object
  return(lgcm)
}

effectiveError <- function(lgcm)
{
  if (!inherits(lgcm,"lgcm")) {
    stop("Not implemented for given class.")
  }
  
  eta <- 1 / ((lgcm$num.timepoints+ lgcm$residual.variance / lgcm$intercept.variance))
  
  rho <-  lgcm$intercept.slope.covariance
  
  zeta <- rho*rho*(lgcm$sumti*lgcm$sumti-lgcm$num.timepoints*lgcm$sumtisq)+
    2*rho*lgcm$sumti*lgcm$residual.variance
  zeta <- zeta / 
    (lgcm$sumtisq*(lgcm$num.timepoints*lgcm$intercept.variance+lgcm$residual.variance)-   
       lgcm$sumti*lgcm$sumti*lgcm$intercept.variance)
 # cat("zeta",zeta)
  lgcm$residual.variance/( lgcm$sumtisq - eta*lgcm$sumti*lgcm$sumti)+ zeta
}

#' Reliability and Precision indices
#' @description A set of indicators for study design efficiency. Effective error, Effective Curve Reliability, and Growth Rate Reliability.
#' @param lgcm A \code{lgcm} object representing a linear latent grwoth curve model
#' @aliases grr effectiveError
ecr <- function(lgcm) 
{
  return (lgcm$slope.variance / (lgcm$slope.variance+effective.error(lgcm)))
}

grr <- function(lgcm) 
{
  err <-  lgcm$residual.variance/( lgcm$sumtisq - (1/lgcm$num.timepoints)*lgcm$sumti*lgcm$sumti)
  return (lgcm$slope.variance / (lgcm$slope.variance+err))
}

