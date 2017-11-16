\name{DS.PostMean}
\alias{DS.PostMean}
\alias{DS.PostMean.reduce}
\alias{EXP.denom}
\alias{EXP.num.1}
\alias{EXP.num.2}
\alias{PP.prt.2}
\alias{weight.fun.beta}
\alias{EXP.score}
\title{
Conditional Expectation for DS Prior
}
\description{ A function that determines the conditional expectation for estimation and posterior means.  
}
\usage{
DS.PostMean(y.i, n.i, g.par, u, LP.par, method = c("LP-Score", "Theta"))
DS.PostMean.reduce(DS.GF.obj)
}
\arguments{
  \item{y.i}{ Integer of number of successes \eqn{y_i} for a single partition
}
  \item{n.i}{ Integer of trial size \eqn{n_i} for a single partition}
  
  \item{g.par}{ Parameters (\eqn{\alpha} and \eqn{\beta}) for parametric beta \eqn{g(\theta)}
}
  \item{u}{ grid: values between 0 and 1 of given length }
  
  \item{LP.par}{vector of LP coefficients}
  
  \item{method}{ \code{LP-Score} sets \eqn{h(\theta_i) = Leg_j(G(\theta_i))}, while \code{Theta} sets \eqn{h(\theta_i) = \theta_i}.}
  
   \item{DS.GF.obj}{object resulting from function \code{DS.prior}}
}
\details{ \code{DS.PostMean.reduce} uses \code{DS.PostMean} to determine the conditional mean for a given set of observations; it requires object generated from function \code{DS.prior}.}

\value{ Returns the LP-Bayes conditional mean estimator for \eqn{(y_i, n_i)} given selected \eqn{h(\theta_i)}. For \code{DS.PostMean.reduce}, returns a vector of the conditional means for all observed values.
}
\references{
Mukhopadhyay, S. and Fletcher, D. (2017), "Bayes Modeling via Goodness of Fit," Technical report.
}
\author{Doug Fletcher
}
\examples{
data(rat)
rat.start <- BetaBinoMLE(rat$y, rat$n)$estimate
rat.ds <- DS.prior(rat, max.m = 4, rat.start)
### Find conditional mean for a specific observation y = 4, n = 14
DS.PostMean(4, 14, rat.ds$g.par, seq(0,1, length.out = 1000), rat.ds$LP.par, method = "Theta")
### Find conditional mean for all observations in a data set
rat.ds.mean <- DS.PostMean.reduce(rat.ds)
}

