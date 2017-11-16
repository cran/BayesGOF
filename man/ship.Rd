\name{ship}
\alias{ship}
\docType{data}
\title{
Portsmouth Navy Shipyard Data
}
\description{
Data represents results of quality-control inspections executed by Portsmouth Naval Shipyard on lots of welding materials.  The data has \eqn{k=5} observations of number of defects \eqn{y} out of the total number of tested \eqn{n=5}.
}
\usage{data("ship")}
\format{
  A data frame of \eqn{(y_i, n_i)} for \eqn{i = 1,...,5}.
  \describe{
    \item{\code{y}}{number of defects found}
    \item{\code{n}}{Total samples tested}
  }
}
\source{ Martz, H. and Lian, M. (1974), "Empirical Bayes estimation of the binomial parameter," Biometrika, 61, 517-523.}
\references{
Sivaganesan, S. and Berger, J. (1993), "Robust Bayesian analysis of the binomial empirical Bayes problem," Canadian Journal of Statistics, 21, 107-119.
}
\keyword{datasets}
