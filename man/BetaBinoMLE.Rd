\name{BetaBinoMLE}
\alias{BetaBinoMLE}
\alias{mom.est}

\title{ Beta-Binomial Parameter Estimation }
\description{
A function that uses maximum likelihood to determine the parametric parameters for \eqn{g(\theta; \alpha, \beta)}.  
}
\usage{
BetaBinoMLE(success, trials, start = NULL, optim.method = "default", 
lower = 0, upper = Inf)
mom.est(success,trials)
}

\arguments{
  \item{success}{Vector containing the number of successes
}
  \item{trials}{Vector containing the total number of trials that correspond to the successes
}
  \item{start}{initial parameters; default is NULL which allows function to determine MoM estimates as initial parameters
}
  \item{optim.method}{optimization method in \code{optim()}{stats}
}
  \item{lower}{ lower bound for parameters; default is 0
}
  \item{upper}{ upper bound for parameters; default is infinity
}
}
\details{
Also available is \code{mom.est}, which returns method of moments estimates for beta parameters without using maximum likelihood.
}
\value{
  \item{estimate}{MLE estimate for beta parameters}
  \item{convergence }{Convergence code from \code{optim()};  0 means convergence}
  \item{loglik}{Loglikelihood that corresponds with MLE estimated parameters}
  \item{initial}{Initial parameters, either user-defined or determined from method of moments}
  \item{hessian}{Estimated Hessian matrix at the given solution}
  }
\references{\url{https://github.com/SupplyFrame/EmpiricalBayesR/blob/master/EmpiricalBayesEstimation.R
}
}

\author{Aleksandar Bradic
}

\examples{
data(rat)
### MLE estimate of alpha and beta
rat.mle <- BetaBinoMLE(rat$y, rat$N)$estimate
rat.mle
### MoM estimate of alpha and beta
rat.mom <- BetaBinoMLE(rat$y, rat$N)$initial
rat.mom
rat.mom.2 <- mom.est(rat$y, rat$N)
rat.mom.2}