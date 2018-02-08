\name{DS.PostMode}
\alias{DS.PostMode}
\alias{DS.mode.map.pgu}
\alias{DS.mode.map.nnu}
\alias{DS.mode.map.bbu}
\alias{DS.PostMode.reduce}
\alias{DS.mode.reduce.pgu}
\alias{DS.mode.reduce.nnu}
\alias{DS.mode.reduce.bbu}

\title{
Mode Finder for DS Prior
}
\description{ A function that determines the posterior mode.  
}
\usage{
DS.PostMode(y.0, n.0 = NULL, g.par, LP.par, B=250, 
            family = c("Normal","Binomial", "Poisson"))
			
DS.PostMode.reduce(DS.GF.obj)
}
\arguments{
  \item{y.0}{ \eqn{y_0} as the number of successes (\code{"Binomial"}), mean (\code{"Normal"}), or count (\code{"Poisson"})
}
  \item{n.0}{ \eqn{n_0} as the trial size for \code{"Binomial"} or \eqn{se.i} for \code{"Normal"}
  }
    \item{g.par}{Vector with estimated parameters for specified conjugate prior distribution \eqn{g} (i.e beta prior: \eqn{\alpha} and \eqn{\beta}; normal prior: \eqn{\mu} and \eqn{\tau^2}; gamma prior:\eqn{\alpha} and \eqn{\beta}). }

  \item{LP.par}{vector of LP coefficients}
  
  \item{B}{ number of desired grid values }
  
  \item{family}{The appropriate conjugate pair for the analysis: \code{"Normal"} is the Normal-Normal, \code{"Binomial"} is the Binomial-Beta pair, and \code{"Poisson"} is the Poisson-Gamma pair.}

  \item{DS.GF.obj}{Object resulting from running DS.prior function on a data set}

}

\value{ Returns the mode for \eqn{\hat{\pi}(\theta_0 | y_0)}. For \code{DS.PostMode.reduce}, returns the mode for \eqn{\hat{\pi}(\theta_i | y_i)}, where \eqn{i = 1,...,k}.
}
\references{
Mukhopadhyay, S. and Fletcher, D. (2018), "Bayesian Modeling via Goodness of Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.

Mukhopadhyay, S. (2017). Large-scale mode identification and data-driven sciences. Electronic Journal of Statistics, 11(1), 215-240.
}
\author{Doug Fletcher, Subhadeep Mukhopadhyay
}
\examples{
data(rat)
rat.start <- gMLE.bb(rat$y, rat$n)$estimate
rat.ds <- DS.prior(rat, max.m = 4, rat.start)
### Find mode for a specific observation y = 4, n = 14
DS.PostMode(4, 14, rat.ds$g.par, rat.ds$LP.par, family = "Binomial")
### Find mode for all observations in a data set
rat.ds.mode <- DS.PostMode.reduce(rat.ds)
}

