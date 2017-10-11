#
# Power simulation replication of Hertzog et al., 2008
# "Evaluating the Power of Latent Growth Curve Models to Detect Individual
# Differences in Change"
# DOI: 10.1080/10705510802338983
#
# see page 549 for parameter values
#
# Statistical power of the generalized variance test (2df) and
# the specific variance test (1df, assuming zero intercept-slope
# correlation) to reject
# a hypothesis of no slope variance.
# We are recreating a line from Figure 5 (p.557) bottom panel of left column 
# (sample size=500 and 4 occasions)
#
# using parallel execution with Snowfall
#
require("OpenMx")
require("sempower")

manifests<-c("x1","x2","x3","x4")
latents<-c("icept","slope")
#errs <- c(1,10,25,100)
err <- 10
icept <- 100
slope <- 25
ti <- c(0,2,4,6) / 20.0

iceptslopefree <- TRUE
N <- 500 # sample size

gcr <- round(icept/(icept+err),2)

cr <- 0
iceptslope <- cr * (sqrt(slope)*sqrt(icept)) # intercept-slope covariance

true.model <- mxModel("LGCM", 
                      type="RAM",
                      manifestVars = manifests,
                      latentVars = latents,
                      mxPath(from="icept",to=c("x1","x2","x3","x4"), free=c(FALSE,FALSE,FALSE,FALSE), 
                             value=c(1.0,1.0,1.0,1.0) , 
                             arrows=1, label=c("icept_TO_x1","icept_TO_x2","icept_TO_x3","icept_TO_x4") ),
                      mxPath(from="slope",to=c("x1", "x2","x3","x4"), free=c(FALSE,FALSE,FALSE,FALSE), 
                             value=c(ti[1],ti[2],ti[3],ti[4]) , arrows=1, 
                             label=c("slope_TO_x1","slope_TO_x2","slope_TO_x3","slope_TO_x4") ),
                      mxPath(from="one",to=c("icept","slope"),free=T,value=c(50,-20)),
                      mxPath(from="icept",to=c("icept"), free=c(TRUE), value=c(icept) ,
                             arrows=2, label=c("sigma_i") ),
                      mxPath(from="icept",to=c("slope"), free=c(iceptslopefree), value=c(iceptslope) ,
                             arrows=2, label=c("sigma_is") ),
                      mxPath(from="slope",to=c("slope"), free=c(TRUE), value=c(slope) , 
                             arrows=2, label=c("sigma_s") ),
                      mxPath(from="x1",to=c("x1"), free=c(TRUE), value=c(err) , arrows=2, label=c("e") ),
                      mxPath(from="x2",to=c("x2"), free=c(TRUE), value=c(err) , arrows=2, label=c("e") ),
                      mxPath(from="x3",to=c("x3"), free=c(TRUE), value=c(err) , arrows=2, label=c("e") ),
                      mxPath(from="x4",to=c("x4"), free=c(TRUE), value=c(err) , arrows=2, label=c("e") )
);

mxModel <- getOpenMxRepresentation(lgcm())

sseq <- seq(-.5,.5,0.1)
result.1df <- rep(NA, length(sseq))
result.2df <- rep(NA, length(sseq))


for (i in 1:length(sseq)) {

  cr <- sseq[i]
  
iceptslope <- cr * (sqrt(slope)*sqrt(icept)) # intercept-slope covariance

temp.true.model <- omxSetParameters(true.model,labels="sigma_is",values=iceptslope)
temp.h0.model <- omxSetParameters(temp.true.model,labels=c("sigma_s","sigma_is"),values=0, free=FALSE)
temp.h1.model <- omxSetParameters(temp.true.model,labels=c("sigma_is"),values=0,free = FALSE)
#simp.1df <- simPowerZeroRestriction(temp.true.model, c("sigma_s"), N=N, 
#                                    repetitions=200, snowfallCpus=6)
simp.1df <- simPower(temp.h0.model,temp.h1.model, temp.true.model,N = N,repetitions = 200)

temp.true.model <- omxSetParameters(true.model,labels="sigma_is",values=iceptslope)
simp.2df <- simPowerZeroRestriction(temp.true.model, c("sigma_s","sigma_is"), N=N, 
                                repetitions=200, snowfallCpus=6)



estpower.2df <- sum(simp.2df$p<0.05)/length(simp.2df$p)
estpower.1df <- sum(simp.1df$p<0.05)/length(simp.1df$p)

result.1df[i] <- estpower.1df
result.2df[i] <- estpower.2df

}


#
# plot diagram
#


lo.1df <- smooth.spline(sseq, result.1df, spar=0.45)
lo.2df <- smooth.spline(sseq, result.2df, spar=0.45)
plot(sseq,result.1df, ylim=c(0,1),type="n",xlab="Intercept/Slope Correlation",ylab="Power")
title("Occasions = 4")
lines(x = sseq, y=result.1df,lty=3)
lines(x = sseq, y=result.2df, lty=3)
lines(x = sseq, y = predict(lo.1df)$y,lw=2)
points(x = sseq, y = predict(lo.1df)$y,pch=1)
lines(x = sseq, y = predict(lo.2df)$y, lw=2)
points(x = sseq, y = predict(lo.2df)$y, pch=4)
abline(h=c(0,0.2,0.4,0.6,0.8,1))
#text(sseq[13], result[13], paste("GCR=",gcr),pos=4)