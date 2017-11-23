#' Estimate Power or Sample Size from a simulation
#'
#' @examples:
#'
#' power(psim, alpha=.05)
#' sampleSize(psim, power=0.8)
#'
#'
#' @param psim A simPower result object.
#' @param alpha Alpha level for power estimation
#' @param power Target power
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

power <- function(psim, alpha=.05, exclude=c())
{
  selector <- rep(TRUE,length(psim$p))
  if ("negstat" %in% exclude) {
	  selector <- psim$statistics>=0
  }

  success <- (psim$p[selector] <alpha)
  
  return( sum(success)/sum(selector))
}
