% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutr.R
\name{kest.gt}
\alias{kest.gt}
\title{Expected number of classes with more than N individuals given parameters \eqn{\theta} and sample size J using the species frequency spectrum}
\usage{
kest.gt(theta, J, N)
}
\arguments{
\item{theta}{parameter \eqn{\theta},  a real value}

\item{J}{is the number of individuals, an integer}

\item{N}{is the minimal abundance, an integer}
}
\value{
k is the expected number of classes \eqn{>N}, a real number
}
\description{
The expected number of species \eqn{>N} is computed using the approximation
based on the species frequency spectrum \eqn{\psi(x)=\theta x^{-1}(1-x)^{\theta-1}}
such that the number of species \eqn{>N} is the integral from N/J to 1 of \eqn{(1-(1-x)^J)\psi(x)}.
This is an approximation of the exact formula, but it turns out to be accurate.
}
\examples{
kest.gt(11.3,234373,50)
}
