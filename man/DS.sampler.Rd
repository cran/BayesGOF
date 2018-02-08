\name{DS.sampler}
\alias{DS.sampler}
\alias{rDS.bbu}
\alias{rDS.pgu}
\alias{rDS.nnu}
\alias{DS.sampler.post}
\alias{rDS.post.nnu}
\alias{rDS.post.pgu}
\alias{rDS.post.bbu}
\title{
Random Data Sampler for DS Prior Distribution
}
\description{
Generates samples of size \eqn{k} from DS\eqn{(G,m)} prior distribution.  
}
\usage{
DS.sampler(k, g.par, LP.par, con.prior)

DS.sampler.post(k, g.par, LP.par, y.0, n.0, con.prior)
}
\arguments{
  \item{k}{Total number of samples requested}
  \item{g.par}{Estimated parameters for specified conjugate prior distribution (i.e beta prior: \eqn{\alpha} and \eqn{\beta}; normal prior: \eqn{\mu} and \eqn{\tau^2}; gamma prior:\eqn{\alpha} and \eqn{\beta}).}
  \item{LP.par}{LP coefficients for DS prior}
  \item{con.prior}{The distribution type of conjugate prior \eqn{g}; either  \code{"Beta"}, \code{"Normal"}, or \code{"Gamma"}.}
  \item{y.0}{Depending on \eqn{g}, \eqn{y_0} is either (i) the sample mean (\code{"Normal"}), (ii) the number of successes (\code{"Beta"}), or (iii) the specific count value (\code{"Gamma"}) for desired posterior distribution(\code{DS.samp.post} only)}
  \item{n.0}{Total number of trials \eqn{n_0} for desired posterior (\code{DS.samp.post} only)}
}

\details{\code{DS.samp.post} uses the same type of sampling to generate random values from a DS posterior distribution.}

\value{
vector of length \eqn{k} containing sampled values from DS prior or DS posterior}
\references{
Mukhopadhyay, S. and Fletcher, D. (2018), "Bayesian Modeling via Goodness of Fit," Technical report, \url{https://arxiv.org/abs/1802.00474 }.

Mukhopadhyay, S. (2017), "Large-Scale Mode Identification and Data-Driven Sciences," Electronic Journal of Statistics, 11, 215-240.

}
\author{Doug Fletcher, Subhadeep Mukhopadhyay}
\examples{
##Extracted parameters from rat.ds object
rat.g.par <- c(2.3, 14.1)
rat.LP.par <- c(0, 0, -0.5)
samps.prior <- DS.sampler(25, rat.g.par, rat.LP.par, con.prior = "Beta")
hist(samps.prior,15)
##Posterior for rat data
samps.post <- DS.sampler.post(25, rat.g.par, rat.LP.par, 
							y.0 = 4, n.0 = 14, con.prior = "Beta")
hist(samps.post, 15)
}