\name{lgcm}
\alias{lgcm}

\title{Create a latent growth curve model}
\description{This function creates a lgcm object representing a linear latent growth curve model}

\usage{
lgcm(timepoints=0, intercept.variance=0, slope.variance=0,
                 residual.variance=0, 
                 intercept.slope.covariance = 0, sample.size = 0)
}

\arguments{
	\item{timepoints}{Vector of time points of measurement occasions}
	\item{sample.size}{Number of cases}
	\item{residual.variance}{Residual variance or inverse of the precision of measurement}
	\item{slope.variance}{Individual differences in slopes}
	\item{intercept.variance}{Individual differences in intercept}
}


\author{Andreas M. Brandmaier}

