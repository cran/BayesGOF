\name{terb}
\alias{terb}
\docType{data}
\title{ Terbinafine trial data
}
\description{ During several studies of an oral antifungal agent, terbinafine, a proportion of the patients in the trial terminated treatment due to some adverse effects.  In the data set, \eqn{y} is the number of terminated treatments and \eqn{n} is the total number of patients in a particular study.
}
\usage{data("terb")}
\format{
 A data frame of \eqn{(y_i, n_i)} for \eqn{i = 1,...,41}.
  \describe{
    \item{\code{y}}{number of patients who terminated treatment early}
    \item{\code{n}}{total number of patients in the clinical trial}
  }
}

\source{ Young-Xu, Y. and Chan, K.A. (2008), "Pooling overdispersed binomial data to estimate event rate," BMC Medical Research Methodology, 8, 58.
}
\keyword{datasets}
