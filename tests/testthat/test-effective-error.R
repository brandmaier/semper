#context("Test effective error")

test_that("error from ECR is identical to analytical solution", {

  library(semper)
  model <- lgcm( timepoints=c(0,1,2,3), intercept.variance = 10,
                 slope.variance = 15, residual.variance = c(10),
                 intercept.slope.covariance = 0)
  
  effective_error_estimate1 <- effectiveError(model)
  
  ecr_estimate <-  ecr(model)
  effective_error_estimate2 <- effective.error.from.ecr( 15, ecr_estimate)
  
  testthat::expect_equal(effective_error_estimate1, effective_error_estimate2)
})
