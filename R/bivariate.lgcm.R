#' Bivariate Latent Growth Curve Model
#' 
#' @param lgcm.x LGCM object for first process (refered to as x)
#' @param lgcm.y LGCM object for second process (refered to as y)
#' @param slope.x.slope.y.covariance Latent covariance of slope of x and slope of y
#' @param icept.x.slope.y.covariance  Latent covariance of intercept of x and slope of y
#' @param icept.y.slope.x.covariance  Latent covariance of slope of x and intercept of y
#' @param icept.x.icept.y.covariance  Latent covariance of intercept of x and intercept of y
#' @param latent.correlations Specify a latent correlation matrix instead of a covariance matrix.
#' @export

bivariateLgcm<-function(lgcm.x, lgcm.y, slope.x.slope.y.covariance=0, icept.x.slope.y.covariance=0,
                         icept.y.slope.x.covariance =0, icept.x.icept.y.covariance=0,
                         latent.correlations=FALSE)
{
result <- list()
result$lgcm.x <- lgcm.x
result$lgcm.y <- lgcm.y
result$slope.x.slope.y.covariance <- slope.x.slope.y.covariance
result$icept.x.slope.y.covariance <- icept.x.slope.y.covariance
result$icept.y.slope.x.covariance <- icept.y.slope.x.covariance
result$icept.x.icept.y.covariance <- icept.x.icept.y.covariance

result$latent.correlations <- latent.correlations
class(result) <- "bivariateLgcm"

return(result)
}