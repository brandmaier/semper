#
# Power simulation replication of Hertzog et al., 2008
# "Evaluating the Power of Latent Growth Curve Models to Detect Individual
# Differences in Change"
# DOI: 10.1080/10705510802338983
#
# see page 549 for parameter values
#
# Statistical power of the generalized variance test (2df) to reject
# a hypothesis of no slope variance.
# We are recreating a line from Figure 1 (p.552)
#
# using parallel execution with Snowfall
#
require("OpenMx")
require("semper")

manifests<-c("x1","x2","x3","x4")
latents<-c("icept","slope")
errs <- c(1,10,25,100)
err <- errs[2]
icept <- 100
slope <- 50
ti <- c(0,2,4,6) / 20.0

iceptslopefree <- TRUE
N <- 500 # sample size e (100,200,500)

mcreps <- 200

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

#sseq <- seq(-1,1,0.1)
#result <- rep(NA, length(sseq))
#
#for (i in 1:length(sseq)) {#
#
#  cr <- sseq[i]
#  
#iceptslope <- cr * (sqrt(slope)*sqrt(icept)) # intercept-slope covariance
#true.model <- omxSetParameters(true.model,labels="sigma_is",values=iceptslope)#
#
#simp <- simPowerZeroRestriction(true.model, c("sigma_s","sigma_is"), N=N, 
#                                repetitions=200, snowfallCpus=7)#
#
#estpower <- sum(simp$p<0.05)/length(simp$p)
#
#result[i] <- estpower
#
#}

simp <- simPowerZeroRestriction(true.model, c("sigma_s","sigma_is"), N=N, 
                                repetitions=mcreps, 
                                steps=list("sigma_is"=seq(-1,1,0.1)*(sqrt(slope)*sqrt(icept)),
                                           "e"=c(1,10,25,100)
                                           ))



# 
# ECR
#
ecrs <- apply(X=simp$combn, 1, function(x) {
  model <- lgcm(timepoints = ti,intercept.variance = icept, slope.variance = slope,
                residual.variance = x[2],intercept.slope.covariance = x[1])
  return(ecr.generic(model,c("slope.variance","intercept.slope.covariance")))
}  )

cor_is <- edata$sigma_is/max(edata$sigma_is)

edata <- cbind(simp$combn, ecrs, cor_is)

pl2 <- ggplot2::ggplot(data=edata, aes(x = cor_is, y=ecrs, group=e)) + 
  geom_line()+
  theme_light() + xlab("Intercept/Slope Correlation")+
  ylab("ECR") + 
  geom_smooth(method = "loess",color="black",se=FALSE)

plot(pl2)

#
#
#

pl <- ggplot2::ggplot(data=edata, aes(x = sigma_is, y=ps, group=e)) + 
  geom_line()+
  theme_light() + xlab("Intercept/Slope Correlation")+
  ylab("Power") + 
 # geom_smooth(method = "loess",color="black",se=FALSE)

plot(pl)

ggsave(filename="hertzog-figure.pdf", pl)


library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
pg <- grid.arrange(pl, pl2, ncol=2)

ggsave(pg, filename="hertzog-figure-table.pdf")

#
# plot diagram
#



title("Occasions = 4")
lo <- smooth.spline(sseq, result, spar=0.35)
plot(sseq,result, ylim=c(0,1),type="l",xlab="Intercept/Slope Correlation",ylab="Power")
lines(x = sseq, y = predict(lo)$y,lw=2)
abline(h=c(0,0.2,0.4,0.6,0.8,1))
text(sseq[13], result[13], paste("GCR=",gcr),pos=4)