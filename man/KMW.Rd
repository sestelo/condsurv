\name{KMW}
\alias{KMW}
\title{Kaplan-Meier weights.}
\description{This function returns a vector with the Kaplan-Meier weights.}
\usage{
KMW(time, status)
}
\arguments{

  \item{time}{ Survival time of the process.}
  \item{status}{Censoring indicator of the survival time of the process; 0 if the survival time is censored and 1 otherwise. }
  
}

\value{ Vector with Kaplan-Meier weights. }

\references{
E. Kaplan and P. Meier. Nonparametric estimation from incomplete observations. Journal of the
American Statistical Association, 53:457-481, 1958.
}

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{PKMW}}}

\examples{
data(colonCS)
obj <- with(colonCS, survCS(time1, event1, Stime, event, age = age))
obj0 <- obj[[1]]
KMW(time = obj0$Stime, status = obj0$event)
}

