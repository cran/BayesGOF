\name{LP.basis.beta}
\alias{LP.basis.beta}
\title{
Determine LP basis functions (beta distributions)
}
\description{Determines the LP basis for a given parametric beta distribution.
}
\usage{
LP.basis.beta(x, g.par, m, ind)
}
\arguments{
  \item{x}{ \code{x} values (integer or vector) from 0 to 1}
  \item{g.par}{ Parameters \eqn{\alpha} and \eqn{\beta} for beta distribution}
  \item{m}{Number of legendre polynomials}
  \item{ind}{Index of desired column (up to m) to return; default is NULL, which returns all m columns}
}
\value{ Matrix with \code{m} columns of values for the LP Score functions evaluated at \code{x}}
\references{
Mukhopadhyay, S. and Fletcher, D. (2017), "Bayes Modeling via Goodness of Fit," Technical report.

Mukhopadhyay,S. and Parzen, E. (2014), "LP Approach to Statistical Modeling," arXiv: 1405.2601.

}
\author{ Doug Fletcher, Subhadeep Mukhopadhyay
}
\keyword{DS Main Functions}