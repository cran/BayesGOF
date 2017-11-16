\name{rDS}
\alias{rDS}
\alias{rDS.posterior}
\title{
Random Data Sampler for DS Prior Distribution
}
\description{
Generates samples of size \eqn{k} from \eqn{DS(G,m)} prior distribution.  
}
\usage{
rDS(k, g.par, LP.par)
rDS.posterior(k, g.par, LP.par, y.i, n.i)
}
\arguments{
  \item{k}{Total number of samples requested}
  \item{g.par}{\eqn{\alpha} and \eqn{\beta} parameters for parametric \eqn{g}}
  \item{LP.par}{LP coefficients for DS prior}
  \item{y.i}{Number of successes \eqn{y_i} for desired posterior (rDS.posterior only)}
  \item{n.i}{Total number of trials \eqn{n_i} for desired posterior (rDS.posterior only)}
}

\details{\code{rDS.posterior} uses same type of sampling to generate random values from a DS posterior distribution.}

\value{
vector of length \eqn{k} containing sampled values from DS prior or DS posterior}
\references{
Mukhopadhyay, S. and Fletcher, D. (2017), "Bayes Modeling via Goodness of Fit," Techincal report.

Mukhopadhyay, S. (2017), "Large-Scale Mode Identification and Data-Driven Sciences," Electronic Journal of Statistics, 11, 215-240.

}
\author{Doug Fletcher, Subhadeep Mukhopadhyay}
\examples{
##Extracted parameters from rat.ds object
rat.g.par <- c(2.3, 14.1)
rat.LP.par <- c(0, 0, -0.5)
samps.prior <- rDS(100, rat.g.par, rat.LP.par)
hist(samps.prior,15)
##Posterior for rat data
samps.post <- rDS.posterior(100, rat.g.par, 
							rat.LP.par, y.i = 4, n.i = 14)
hist(samps.post, 15)
}