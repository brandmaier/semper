#
# Determine the power of a test of mean differences in two groups 
# of equal sample size
#
require("OpenMx")
require("semper")

# sample size
N <- 50
# difference in means
mu.delta <- .3

manifests<-c("x1","x2")
latents<-c()
model <- mxModel("True Model", 
                 type="RAM",
                 manifestVars = manifests,
                 latentVars = latents,
                 mxPath(from="one",to=c("x1","x2"), free=c(TRUE,TRUE), value=c(0,mu.delta) , arrows=1, label=c("mu_x1","mu_x2") ),
                 mxPath(from="x1",to=c("x1"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_x1") ),
                 mxPath(from="x2",to=c("x2"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_x2") )
);

model.alt <- mxModel("H0 Model", 
                 type="RAM",
                 manifestVars = manifests,
                 latentVars = latents,
                 mxPath(from="one",to=c("x1","x2"), free=c(TRUE,TRUE), value=c(1.0,1.0) ,
                        arrows=1, label=c("mu","mu") ),
                 mxPath(from="x1",to=c("x1"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_x1") ),
                 mxPath(from="x2",to=c("x2"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_x2") )
);


# generate some data

x1 <- rnorm(n = N, mean=0, sd=1)
x2 <- rnorm(n = N, mean=mu.delta, sd=1)

ttest.result <- power.t.test(n=N,delta=mu.delta,sd=1)

sp.result <- simPower(model.alt, model, model, N=N, repetitions=200)

print(ttest.result)
print(sp.result)
