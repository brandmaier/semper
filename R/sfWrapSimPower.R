

sfWrapSimPower <- function(repetitions,h0Model,h1Model, populationModel, 
                           N.range, keep.models)
{
  return(simPower(h0Model,h1Model, populationModel,N.range,repetitions,keep.models ))
}
