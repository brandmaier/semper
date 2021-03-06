#' Estimate Power or Sample Size from a simulation
#'
#' @examples 
#'
#' power(psim, alpha=.05)
#' sampleSize(psim, power=0.8)
#'
#'
#' @param psim A simPower result object.
#' @param alpha Alpha level for power estimation
#' @param power Target power
#' @exclude Vector of exclusion criteria. If "negstat" then negative values of statistics are excluded.
#'
#' @return Return either the sample size if power was set with \code{estimatePower} or the power if sample size was fixed with \code{estimateSampleSize}.
#'
#' @author Andreas M. Brandmaier
#'
#' @aliases sampleSize
#'
#' @seealso \code{\link{simPower}}
#'
#' @export


power<-function (psim, alpha = 0.05, exclude = c("na")) 
{
  if ("na" %in% exclude) {
    psim$p <- psim$p[!is.na(psim$p)]
    psim$z <- psim$z[!is.na(psim$p)]
    psim$estimates <- psim$estimates[!is.na(psim$p)]
    psim$N <- psim$N[!is.na(psim$p)]
  } else {
    psim$p[is.na(psim$p)]<-FALSE
  }
  
  selector <- rep(TRUE, length(psim$p))
  if ("negstat" %in% exclude) {
    selector <- psim$statistics >= 0
  }
  success <- (psim$p[selector] < alpha)
  return(sum(success)/sum(selector))
}