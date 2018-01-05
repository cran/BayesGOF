\name{DS.prior}
\alias{DS.prior}
\alias{plot.DS_GF}
\alias{print.DS_GF}
\alias{Reduce.LP.coef}
\alias{DS.GF.obj}
\alias{LP.smooth}
\title{DS Prior: Diagnostics and Estimation
}
\description{A function to estimate DS\eqn{(G,m)} prior model along with diagnostics.
}
\usage{
DS.prior(yn.df, max.m = 8, start.par, iter.c = 200,
			B = 1000, smooth.crit = "BIC")
}
\arguments{
  \item{yn.df}{dataframe that contains the \eqn{k} pairs of successes \eqn{y} and the corresponding total number of trials \eqn{n}}
  \item{max.m}{Integer value that reflects the maximum amount of flexibility to deviate from \eqn{g}; default value is set at 8 (appropriate for most data sets)}
  \item{start.par}{Vector containing the parameters \eqn{\alpha} and \eqn{\beta} for parametric beta distribution \eqn{g}}
  \item{iter.c}{Integer value that gives the maximum number of iterations allowed for convergence; default is 200}
  \item{B}{Integer value for number of grid points used for distribution output; default is 1000}
  \item{smooth.crit}{User selects either \code{"BIC"} or \code{"AIC"} as criteria to both determine optimal \eqn{m} and smooth final LP parameters; default is \code{"BIC"}}
}
\details{ The dataframe and starting parameters for parametric beta distribution \eqn{g} are always required.  Function can take \eqn{m=0} and will return the Empirical Bayes estimate with given starting parameters for the result.  Returns an object of class \code{DS.GF.obj}; this object can be used with plot command to plot the U-function (\code{Ufunc}), Deviance Plots (\code{mDev}), and DS-G comparison (\code{DS_G}).
}
\value{
  \item{LP.par}{\eqn{m} smoothed LP-Fourier coefficients, where \eqn{m} is determined by maximum deviance}
  \item{g.par}{Parametric parameters for \eqn{g(\theta; \alpha,\beta)}; same values as given by user for start.par}
  \item{LP.max.smt}{All potential smoothed LP-Fourier coefficients, where the number is the same as \code{max.m}}
  \item{LP.max.uns}{All potential LP-Fourier coefficients prior to smoothing, where the number is the same as \code{max.m}}
  \item{m.val}{The value of \eqn{m} (less than or equal to the maximum \eqn{m} from user) that has the maximum deviance and represents the optimal number of LP-Fourier coefficients}
  \item{obs.data}{Dataframe for original observed data user provided as \code{yn.df}}
  \item{prior.data}{Dataframe that contains B values for \eqn{\theta} between 0 and 1 (\code{$theta.vals}), the corresponding density values of the parametric prior \eqn{g} (\code{$parm.prior}), and the corresponding density values of \eqn{\pi(\theta)} (\code{$ds.prior})}
  \item{UF.data}{Dataframe that contains values required for plotting the U-function}
  \item{dev.df}{Dataframe that contains deviance values for values of \eqn{m} up to \code{max.m}}
  \item{sm.crit}{Smoothing criteria; either \code{"BIC"} or \code{"AIC"}}
  }
\references{
Mukhopadhyay, S. and Fletcher, D. (2017), "Bayes Modeling via Goodness of Fit," Technical report.

Mukhopadhyay, S. (2017), "Large-Scale Mode Identification and Data-Driven Sciences," Electronic Journal of Statistics, 11, 215-240.

}
\author{Doug Fletcher, Subhadeep Mukhopadhyay}
\examples{
data(rat)
rat.start <- BetaBinoMLE(rat$y, rat$n)$estimate
rat.ds <- DS.prior(rat, max.m = 4, rat.start)
rat.ds
plot(rat.ds, plot.type = "Ufunc")
plot(rat.ds, plot.type = "DSg")
plot(rat.ds, plot.type = "mDev")
}
\keyword{DS Main Functions}