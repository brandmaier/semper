\name{estimatePower}
\alias{estimatePower}
\alias{estimateSampleSize}

\title{Estimate Power or Sample Size from a simulation}

\usage{

estimatePower(psim, alpha=.05)
estimateSampleSize(psim, power=0.8)

}

\arguments{
\item{psim}{A simPower result object.}

\item{alpha}{Alpha level for power estimation}
\item{power}{Target power}
}

\value{
Return either the sample size if power was set with \code{estimatePower} or the power if sample size was fixed with \code{estimateSampleSize}.

}

\author{Andreas M. Brandmaier}