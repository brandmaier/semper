# install and load package
#devtools::install_github("brandmaier/semper")
library(semper)

#
# This is the model specification. To obtain 1df-Wald-type reliability
# and effective error, set intercept.variance to something very large (10000) and intercept.slope.covariance=0
#
model <- lgcm( timepoints=c(0,1,2,3), intercept.variance = 10,
      slope.variance = 15, residual.variance = c(10,10,20,30),
      intercept.slope.covariance = 0)

# estimate ECR and reliability
ecr.estimate <- ecr.generic(model)
eff2d <- effective.error.from.ecr(model$slope.variance, ecr.estimate)

# output results
cat("ECR: ", ecr.estimate,"\n")
cat("Effective error (LR;2df): ", eff2d,"\n")


