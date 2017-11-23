#' @export

simPowerSpecificVarianceTest <- function(model, repetitions=500, sample.size=NULL, full.output=FALSE) {

 trueModel <- getOpenMxRepresentation(model)

 h1Model <- trueModel

 h1Model <- omxSetParameters(h1Model, labels="interceptslopecovariance",free = FALSE,values = 0)
 h0Model <- h1Model
 h0Model <- omxSetParameters(h0Model, labels="slopevariance",free = FALSE,values = 0)


 if (is.null(sample.size)) {
  N <- model$sample.size
 } else {
  N <- sample.size
 }
 

 simulation <- simPower(trueModel, h1Model=h1Model, 
h0Model=h0Model,N = N, repetitions = repetitions)

 if (full.output) {
return(simulation)
	} else {
 return(power(simulation))
}

}