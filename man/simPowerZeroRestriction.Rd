\name{simPowerZeroRestriction}
\alias{simPowerZeroRestriction}
\title{Simulate Power for a null hypothesis involving one or more zero-constraints.}

\description{Monte-Carlo-simulate power for a given model to reject a null hypothesis involving one or more zero-constraints.}

\usage{
simPowerZeroRestriction(trueModel, restrictions, N, repetitions, keepModels=F, 
                                  snowfallCpus=-1)
}

\arguments{
	\item{trueModel}{An OpenMx model representing the population model.}
	\item{restrictions}{A list of parameters to be restricted to zero under H0.}
	\item{N}{Sample size}
	\item{repetitions}{The number of Monte-Carlo-trials.}
	\item{keepModels}{Store the complete intermediate results in form of the fitted models. Only recommended for debugging purposes. }
	\item{snowfallCpus}{ By default (-1), the simulation is executed sequentally. For any value larger than zero, the execution is executed in parallel using \code{snowfall} package with the specified number of cpus.}
}

\author{Andreas M. Brandmaier}

\details{
	This is a wrapper function to the generic \code{simpower} function to obtain power for a given model to reject a null hypothesis involving one or more zero-constraints by a Monte Carlo simulation.
	 

}
