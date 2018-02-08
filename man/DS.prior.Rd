\name{DS.prior}
\alias{DS.prior}
\alias{DS.prior.nnu}
\alias{DS.prior.pgu}
\alias{DS.prior.bbu}
\alias{plot.DS_GF}
\alias{print.DS_GF}
\alias{Reduce.LP.coef.univ}
\alias{DS.GF.obj}
\alias{LP.smooth}
\title{Prior Diagnostics and Estimation
}
\description{A function that generates the diagnostic uncertainty function (\code{Ufunction}) and estimates DS\eqn{(G,m)} prior model.
}
\usage{
DS.prior(input, max.m = 8, g.par, iters = 200, 
         family = c("Normal","Binomial", "Poisson"),
         B = 1000, smooth.crit = "BIC")
}
\arguments{
  \item{input}{For \code{"Binomial"}, a dataframe that contains the \eqn{k} pairs of successes \eqn{y} and the corresponding total number of trials \eqn{n}.  For \code{"Normal"}, a dataframe that has the \eqn{k} means \eqn{y_i} in the first column and their respective standard errors \eqn{s_i} in the second.  For the \code{"Poisson"}, a vector of that includes the untabled count data.}
  \item{max.m}{Integer value that reflects the maximum amount of flexibility to deviate from \eqn{g}; default value is set at 8 (appropriate for most data sets)}
  \item{g.par}{Vector with estimated parameters for specified conjugate prior distribution \eqn{g} (i.e beta prior: \eqn{\alpha} and \eqn{\beta}; normal prior: \eqn{\mu} and \eqn{\tau^2}; gamma prior: \eqn{\alpha} and \eqn{\beta}). }
  \item{iters}{Integer value that gives the maximum number of iterations allowed for convergence; default is 200}
  \item{family}{The appropriate conjugate pair for the analysis: \code{"Normal"} is the Normal-Normal, \code{"Binomial"} is the Binomial-Beta pair, and \code{"Poisson"} is the Poisson-Gamma pair.}
  \item{B}{Integer value for number of grid points used for distribution output; default is 1000}
  \item{smooth.crit}{User selects either \code{"BIC"} or \code{"AIC"} as criteria to both determine optimal \eqn{m} and smooth final LP parameters; default is \code{"BIC"}}
}
\details{ Function can take \eqn{m=0} and will return the parametric empirical Bayes estimate with given starting parameters for the result.  Returns an object of class \code{DS.GF.obj}; this object can be used with plot command to plot the U-function (\code{Ufunc}), Deviance Plots (\code{mDev}), and DS-G comparison (\code{DS_G}).
}
\value{
  \item{LP.par}{\eqn{m} smoothed LP-Fourier coefficients, where \eqn{m} is determined by maximum deviance}
  \item{g.par}{Parameters for \eqn{g}}
  \item{LP.max.smt}{All smoothed LP-Fourier coefficients, where the number is the same as \code{max.m}}
  \item{LP.max.uns}{All  LP-Fourier coefficients prior to smoothing, where the number is the same as \code{max.m}}
  \item{m.val}{The value of \eqn{m} (less than or equal to the maximum \eqn{m} from user) that has the maximum deviance and represents the appropriate number of LP-Fourier coefficients}
  \item{obs.data}{Observed data provided by user for \code{input}}
  \item{prior.fit}{Fitted values for the estimated prior}
  \item{UF.data}{Dataframe that contains values required for plotting the U-function}
  \item{dev.df}{Dataframe that contains deviance values for values of \eqn{m} up to \code{max.m}}
  \item{sm.crit}{Smoothing criteria; either \code{"BIC"} or \code{"AIC"}}
  \item{fam}{The user-selected family.}
  }
\references{
Mukhopadhyay, S. and Fletcher, D. (2018), (2018), "Bayesian Modeling via Goodness of Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.

Mukhopadhyay, S. (2017), "Large-Scale Mode Identification and Data-Driven Sciences," Electronic Journal of Statistics, 11, 215-240.

}
\author{Doug Fletcher, Subhadeep Mukhopadhyay}
\examples{
data(rat)
rat.start <- gMLE.bb(rat$y, rat$n)$estimate
rat.ds <- DS.prior(rat, max.m = 4, rat.start, family = "Binomial")
rat.ds
plot(rat.ds, plot.type = "Ufunc")
plot(rat.ds, plot.type = "DSg")
plot(rat.ds, plot.type = "mDev")
}
\keyword{DS Main Functions}