\name{DS.PostMean}
\alias{DS.PostMean}
\alias{DS.PostMean.pge}
\alias{DS.PostMean.pgu}
\alias{DS.PostMean.nnu}
\alias{DS.PostMean.bbu}
\alias{DS.PostMean.reduce}
\alias{DS.PostMean.reduce.pgu}
\alias{DS.PostMean.reduce.nnu}
\alias{DS.PostMean.reduce.bbu}
\alias{EXP.denom}
\alias{EXP.num.1}
\alias{EXP.num.2}
\alias{ConMean.prt.2.pg}
\alias{ConMean.prt.2.nn}
\alias{ConMean.prt.2.bb}
\alias{weight.fun.univ}
\alias{EXP.score}
\alias{lambda.i}
\title{
Posterior Expectation for DS Prior
}
\description{ A function that determines the posterior expectation \eqn{E(\theta_0 | y_0)} for new study \eqn{y_0}.  
}
\usage{
DS.PostMean(y.0, n.0 = NULL, g.par, LP.par, B =250, 
            family = c("Normal","Binomial", "Poisson"), e.0= NULL)
			
DS.PostMean.reduce(DS.GF.obj)
}
\arguments{
  \item{y.0}{ The new study's specific \eqn{y_0};  \eqn{y_0} is the number of successes for (\code{"Binomial"}), the mean for (\code{"Normal"}), or count for (\code{"Poisson"})
}
  \item{n.0}{  For \code{"Binomial"}, \code{n.0} is the trial size for the new study; For \code{"Normal"}, \code{n.0} is the standard error for the new study \code{y.0}   }
  
  \item{e.0}{ The study-specific exposure in the case of \code{"Poisson"} with exposure.}
  \item{g.par}{Vector with estimated parameters for specified conjugate prior distribution \eqn{g} (i.e beta prior: \eqn{\alpha} and \eqn{\beta}; normal prior: \eqn{\mu} and \eqn{\tau^2}; gamma prior: \eqn{\alpha} and \eqn{\beta}). }

   \item{B}{ number of desired grid values }
  
  \item{LP.par}{vector of LP coefficients}
  
   \item{family}{The appropriate conjugate pair for the analysis: \code{"Normal"} is the Normal-Normal, \code{"Binomial"} is the Binomial-Beta pair, and \code{"Poisson"} is the Poisson-Gamma pair.
   }
   \item{DS.GF.obj}{Object resulting from running DS.prior function on a data set}

}

\value{ \code{DS.PostMean} returns \eqn{E(\theta_0 | y_0)}, while \code{DS.PostMean.reduce} returns \eqn{E(\theta_i | y_i)} for \eqn{i = 1,..., k}.
}
\references{
Mukhopadhyay, S. and Fletcher, D. (2018), "Bayesian Modeling via Goodness of Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.
}
\author{Doug Fletcher
}
\examples{
data(rat)
rat.start <- gMLE.bb(rat$y, rat$n)$estimate
rat.ds <- DS.prior(rat, max.m = 4, rat.start, family = "Binomial")
### Find conditional mean for a specific observation y = 4, n = 14
DS.PostMean(4, 14, rat.ds$g.par, rat.ds$LP.par, family = "Binomial")
### Find conditional mean for all observations in a data set
rat.ds.mean <- DS.PostMean.reduce(rat.ds)
}

