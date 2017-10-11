icc2 <- function(model) {
  model$intercept.variance/(model$intercept.variance+model$residual.variance/model$num.timepoints)
}

