\name{gLP.basis}
\alias{gLP.basis}
\alias{LP.basis.beta}
\alias{LP.basis.norm}
\alias{LP.basis.gamma}
\title{
Determine LP basis functions for conjugate prior distribution \eqn{g}
}
\description{Determines the LP basis for a given parametric conjugate prior distribution.
}
\usage{
gLP.basis(x, g.par, m, ind, con.prior)
}
\arguments{
  \item{x}{ \code{x} values (integer or vector) from 0 to 1}
  \item{g.par}{ Estimated parameters for specified conjugate prior distribution (i.e beta prior: \eqn{\alpha} and \eqn{\beta}; normal prior: \eqn{\mu} and \eqn{\tau^2}; gamma prior:\eqn{\alpha} and \eqn{\beta}).}
  \item{m}{Number of legendre polynomials}
  \item{ind}{Index of desired column (up to m) to return; default is NULL, which returns all m columns}
  \item{con.prior}{Specified conjugate prior distribution for basis functions.  Options are \code{"Beta"}, \code{"Normal"}, and \code{"Gamma"}}
}
\value{ Matrix with \code{m} columns of values for the LP Score functions evaluated at \code{x}}
\references{
Mukhopadhyay, S. and Fletcher, D. (2018), "Bayesian Modeling via Goodness of Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.

Mukhopadhyay,S. and Parzen, E. (2014), "LP Approach to Statistical Modeling," arXiv: 1405.2601.

}
\author{ Doug Fletcher, Subhadeep Mukhopadhyay
}
\keyword{DS Main Functions}