\name{simPower}
\alias{simPower}
\alias{print.simPower}
\alias{plot.simPower}

\title{Monte-Carlo simulate power}

\description{
Monte-Carlo-simulate power for a given model to reject a null hypothesis of a set of restrictions.
}

\details{
 Given a true model, data is repeatedely simulated
 and the power to detect this difference is calculated for a
 chosen sample size. Alternatively, the required sample size can be
 calculated for a chosen power.
 A simplified wrapper for the simulation of zero constraints is \code{simPowerZeroRestriction}.

}

\usage{

simPower(h0Model,h1Model, populationModel, N=c(20,100),repetitions=100, keepModels=F)


plot(psim, prng=seq(0.4,0.95,0.05), target.power=NULL, lw=2, lty=1,
                           xlab="sample size", ylab="statistical power",
                           main="Monte Carlo Power Simulation", add=F,...)
						   
print(psim, ...)						  

}

\arguments{
\item{h0Model}{Model representing the null hypothesis}
\item{h1Model}{Model representing the alternative hypothesis}
\item{populationModel}{True population model used to generate data from.}
\item{N}{Sample size. For RAM models, this is either (1) a scalar representing a fixed sample size, (2) a list with two items representing a minimum and a maximum for randomly drawn sample sizes (when sample size is the dependent variable), or (3) for multiple group models, this is a named list of sample sizes for the sub models where the names must match the submodels' names.}

\item{psim}{SimPower result object}
\item{prng}{Range for plotting p values}
\item{target.power}{If not null, draw a line to indicate sample size for a specified target power in a power curve plot.}
\item{lw}{line width for power curve}
\item{lty}{line type for power curve}
˜
}

\author{Andreas M. Brandmaier}