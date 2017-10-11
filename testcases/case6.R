require("sempower")
#
# test
model <- sempower::lgcm(timepoints=c(0,1,2,3),intercept.variance = 1,
                        residual.variance = 10,slope.variance = 5,intercept.slope.covariance = 0,
                        sample.size = 100)


ecr.specific <- sempower::ecr(model)

ecr.general <- sempower::ecr.generic(model,"slope.variance")


stopifnot( abs(ecr.specific - ecr.general) <1e-5 )

eff.err.general <- sempower::effective.error.from.ecr(ecr.general,abs.effect = model$slope.variance)

eff.err.specific <- effective.error(model)

stopifnot((eff.err.general-eff.err.specific)<1e-5)