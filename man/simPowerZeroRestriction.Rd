% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simPowerZeroRestriction.R
\name{simPowerZeroRestriction}
\alias{simPowerZeroRestriction}
\title{Simulate Power for a null hypothesis involving one or more zero-constraints.}
\usage{
simPowerZeroRestriction(trueModel, restrictions, N, repetitions,
  keepModels = F, steps = NULL, cluster = NULL)
}
\arguments{
\item{trueModel}{An OpenMx model representing the population model.}

\item{restrictions}{A list of parameters to be restricted to zero under H0.}

\item{N}{Sample size}

\item{repetitions}{The number of Monte-Carlo-trials.}

\item{keepModels}{Store the complete intermediate results in form of the fitted models. Only recommended for debugging purposes.}

\item{steps}{A list of parameter names and values for simulating across multiple conditions}

\item{cluster}{A cluster object for parallel execution}
}
\description{
Monte-Carlo-simulate power for a given model to reject a null hypothesis involving one or more zero-constraints.
}
\details{
This is a wrapper function to the generic \code{simpower} function to obtain power for a given model to reject a null hypothesis involving one or more zero-constraints by a Monte Carlo simulation.
}

