\name{DS.micro.inf}
\alias{DS.micro.inf}
\alias{plot.DS_GF_micro}
\alias{print.DS_GF_micro}
\alias{DS.GF.micro.obj}
\alias{DS.mode.map}
\alias{DS.mode.reduce}
\title{ MicroInference for DS Prior Objects}
\description{ Provide a nonparametric elastic Bayes estimate for a specific \eqn{\theta_i} given a new study, \eqn{(y_i, n_i)}.
}
\usage{
DS.micro.inf(DS.GF.obj, y.i, n.i)

}
\arguments{
  \item{DS.GF.obj}{Object resulting from running DS.prior function on a data set}
  \item{y.i}{ Number of success \eqn{y_i} for new study}
  \item{n.i}{ Trial size \eqn{n_i} for new study}
  }
\details{Returns an object of class \code{DS.GF.micro} that can be used in conjunction with plot command to plot the DS posterior distribution for the new study.}
\value{
  \item{post.data}{Dataframe that contains values for \eqn{\theta}, their corresponding density values for the parametric posterior density, and the corresponding density values for the \eqn{\pi(\theta)} posterior density.}
  \item{post.mean}{Posterior mean for \eqn{\pi(\theta_0 | y_0,n_0)}}
  \item{post.mode}{Posterior mode for \eqn{\pi(\theta| y_0,n_0)}}
  \item{study}{User-provided \eqn{y_0} and \eqn{n_0}}
}
\references{
Mukhopadhyay, S. and Fletcher, D. (2017), "Bayes Modeling via Goodness of Fit," Technical report.
}
\author{Doug Fletcher, Subhadeep Mukhopadhyay}
\examples{
### MicroInference for Naval Shipyard Data: sample where y = 0 and n = 5
data(ship)
ship.ds <- DS.prior(ship, max.m = 2, c(.5,.5))
ship.ds.micro <- DS.micro.inf(ship.ds, y.i = 0, n.i = 5)
ship.ds.micro
plot(ship.ds.micro)
}
\keyword{DS Main Functions}
