\name{bladder4CS}
\alias{bladder4CS}
\docType{data}
\title{Bladder Cancer Recurrences.}
\description{
  bladder4CS is a data frame with 9 variables and 85 observations.
  Data on recurrences of bladder cancer, used by many people to demonstrate methodology for recurrent event modelling.}

\usage{data("bladder4CS")}

\format{
  A data frame with 85 observations on the following 9 variables.
  Below a brief description is given for some of these variables.
  \describe{
    \item{t1}{Time to recurrence/censoring, whichever occurs first.}
    \item{e1}{Recurrence/censoring indicator (recurrence=1, alive=0) for the first time (t1).}
    \item{t2}{Time to censoring/death, whichever occurs first.}
    \item{e2}{Death/censoring indicator (death=1, alive=0) for the second time (t2)}
    \item{t3}{Time to recurrence/censoring, whichever occurs first.}
    \item{e3}{Death/censoring indicator (death=1, alive=0) for the third time (t3)}
    \item{t4}{Time to recurrence/censoring, whichever occurs first.}
    \item{e4}{Death/censoring indicator (death=1, alive=0) for the fourth time (t4)}
    \item{size}{Size (cm) of largest initial tumour.}
    }
}



\references{

  Andrews, D.F., Hertzberg, A.M. (1985), DATA: A Collection of Problems from Many
  Fields for the Student and Research Worker, New York: Springer-Verlag.

  Wei, L.J., Lin, D.Y., Weissfeld, L. (1989), Regression analysis of
  multivariate incomplete failure time data by modeling marginal
  distributions. Journal of the American Statistical Association, 84.
}

\examples{
  data(bladder4CS)
  head(bladder4CS)
}

\keyword{datasets}
