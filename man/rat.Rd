\name{rat}
\alias{rat}
\docType{data}
\title{ Rat Tumor Data
}
\description{Incidence of endometrial stromal polyps in \eqn{k=70} studys of female rats in control group of a 1977 study on the carcinogenic effects of a diabetic drug phenformin.  For each of the \eqn{k} groups, \eqn{y} represents the number of rats who developed the tumors out of \eqn{n} total rats in the group.
}
\usage{data("rat")}
\format{
  A data frame of \eqn{(y_i, n_i)} for \eqn{i = 1,...,70}.
  \describe{
    \item{\code{y}}{number of female rats in a study who developed polyps/tumors}
    \item{\code{n}}{total number of rats in the particular study}
  }
}

\source{ 
National Cancer Institute (1977), "Bioassay of phenformin for possible carcinogenicity," Technical Report No. 7.
}
\references{
Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A., and Rubin, D.B. (2013), Baysian Data Analysis, Third Edition, Chapman & Hall/ CRC Texts in Statistical Science, Taylor & Francis.
Tarone, R.E. (1982), "The use of historical control information in testing for a trend in proportions," Biometrics, 215-220.
}
\keyword{datasets}
