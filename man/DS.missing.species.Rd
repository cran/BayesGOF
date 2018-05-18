\name{DS.missing.species}
\alias{DS.missing.species}

\title{DS Missing Species}
\description{Provides an expected number of new distinct observations, given a ratio of \eqn{m} additional samples to the \eqn{n} samples already observed.
}
\usage{
DS.missing.species(DS.GF.obj, t.ratio, find.se = FALSE, se.reps = 50)

}
\arguments{
  \item{DS.GF.obj}{Object resulting from running DS.prior function on a data set; applies only to those of the Poisson family.}
  \item{t.ratio}{The extrapolation factor: \eqn{\frac{m}{n}}.}
  \item{find.se}{If \code{TRUE}, will use bootstrap to find the standard error of the estimate; default is \code{FALSE}.}
  \item{se.reps}{The number of bootstrap iterations used in finding the standard error.}
  }

\details{Function utilizes DS priors in the following equation:\deqn{E(t)=N \int_0^{\infty} e^{-\theta} \frac{1-e^{-\theta t}}{1 - e^{-\theta}} \pi(\theta) d \theta.}  For more details on this equation, see Efron (2017).}

\value{
  \item{ms.val}{Expected number of new distinct observations in the user-provided \code{t.ratio}.}
  \item{boot.se}{The bootstrap standard error of the point estimate.}
  \item{boot.ms.val}{Vector of bootstrapped point estimates used in finding the standard error.}
  \item{t.ratio}{The user-provide ratio for the calculated expectation.}
  
  }
\references{
Efron, B., 2017. "Bayes, Oracle Bayes, and Empirical Bayes," \url{http://statweb.stanford.edu/~ckirby/brad/papers/2017BayesOBayesEBayes.pdf }.

Mukhopadhyay, S. and Fletcher, D., 2018. "Bayesian Modeling via Goodness-of-Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.
}
\author{Doug Fletcher}
\examples{
data(CorbBfly)
bfly.ds <- DS.prior(CorbBfly, max.m = 7, c(0.0899, 104.395), family = "Poisson")
bfly.ds.ms <- DS.missing.species(bfly.ds, t.ratio = .5)
bfly.ds.ms$ms.val
}

