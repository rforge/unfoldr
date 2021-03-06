% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spheroid.R
\name{binning1d}
\alias{binning1d}
\title{Binning numeric values}
\usage{
binning1d(x, bin, na.rm = FALSE)
}
\arguments{
\item{x}{numeric values to be binned}

\item{bin}{non-decreasingly sorted breaks vector}

\item{na.rm}{logical, default \code{FALSE}, whether to remove missing values, including \code{NaN} in \code{x}}
}
\value{
Vector of count data
}
\description{
Vector of count data
}
\details{
The function provides basic binning (grouping) of numeric values into
classes defined by the breaks vector \code{bin}. The values are binned
according to \eqn{bin[i[j]]<x[j]\leq bin[i[j]+1]} for intervals \eqn{i=1,...,N-1}
and \code{length(bin)=N} of values \eqn{x[j]}, \eqn{j=1,...,|x|}. If \eqn{x[j] > bin[N]} or \eqn{x[j] < bin[1]} then \eqn{x[j]}
is not counted at all.
}
\examples{
	x <- runif(100,0,1)
	bin <- seq(0,1,by=0.1)
	binning1d(x,bin)

}
\author{
M. Baaske
}
