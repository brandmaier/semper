
gcr <- function(model, i=1)
{
if (!inherits(model,"lgcm"))
{
stop("Error!")
}

if (i<1 | i > length(model$timepoints)) {
  return(NA)
}
  
loading <- model$timepoints[i]
  
effect <- model$intercept.variance + loading^2*model$slope.variance+2*loading*model$intercept.slope.covariance
error <- model$residual.variance

return ( effect/(effect+error))

}