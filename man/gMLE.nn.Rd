\name{gMLE.nn}
\alias{gMLE.nn}

\title{ Normal-Normal Parameter Estimation }
\description{
A function that uses maximum likelihood, REML, or method of moments to determine the parametric parameters for \eqn{g\sim}Normal\eqn{(\mu, \tau^2)}.  
}
\usage{
gMLE.nn(value, se, fixed = FALSE, method = c("DL","REML","MoM"))
}

\arguments{
  \item{value}{vector of values
}
  \item{se}{standard error for each value
}
  \item{fixed}{When \code{FALSE}, treats the input as if from a random effects model; otherwise, will treat it as if it a fixed effect.}
  \item{method}{Determines the method to find \eqn{\tau^2}: \code{"DL"} uses Dersimonian and Lard technique, \code{"REML"} uses restricted maximum likelihood, and \code{"MoM"} uses a method of moments technique.
}
}

\value{
  \item{estimate}{vector of estimates for Normal prior parameters.}
  \item{mu.hat}{Estimate for \eqn{\hat{\mu}}.}
  \item{tau.sq}{Estimate for \eqn{\hat{\tau}^2}.}
  \item{method}{User-selected method.}
  }
\references{
Marin-Martinez, F. and Sanchez-Meca, J. (2010), "Weighting by inverse variance or by sample size in random-effects meta-analysis," Educational and Psychological Measurement, 70(1), 56-73.

Brown, L.D. (2008), "In-season prediction of batting averages: A field test of empirical Bayes and Bayes methodologies," The Annals of Applied Statistics, 113-152.
}

\author{Doug Fletcher
}

\examples{
data(ulcer)
### MLE estimate of alpha and beta
ulcer.mle <- gMLE.nn(ulcer$y, ulcer$se, method = "DL")$estimate
ulcer.mle
ulcer.reml <- gMLE.nn(ulcer$y, ulcer$se, method = "REML")$estimate
ulcer.reml
}