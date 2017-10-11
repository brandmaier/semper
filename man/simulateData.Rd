\name{simulateData}
\alias{simulateData}

\title{Simulate Data}
\description{Simulating data from a SEM by drawing from the multivariate normal as described by the model covariance matrix and means vector.}

\usage{
	simulateData(model, N)
}

\arguments{
	\item{model}{An OpenMx model.}
	\item{N}{Desired sample size.}
}

\value{A \code{data.frame} with columns relating to the observed variables of the \code{model} and rows to simulated observations. }

\author{Andreas M. Brandmaier}

