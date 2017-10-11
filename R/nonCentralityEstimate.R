nonCentralityEstimate <- function(model.full, model.restricted,N)
{
Sigma2 <- modelPredictedCov(model.full);
Sigma1 <- modelPredictedCov(model.restricted);
numobs <- dim(Sigma1)[1]
lr <- N*(log(det(Sigma1))
-log(det(Sigma2))+sum(diag(solve(Sigma1)%*%Sigma2)) - numobs)

return(lr);
}