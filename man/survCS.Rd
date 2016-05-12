\name{survCS}
\alias{survCS}
\title{Create a survCS object.}
\description{Creates a "survCS" object, usually used as input to other functions.}
\usage{
survCS(time1, event1, Stime, event, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time1}{First time or censoring time.}
  \item{event1}{Indicator of the first time; 0 if the first time is censored and 1 otherwise.}
  \item{Stime}{The total time of the process.}
  \item{event}{Censoring indicator of the survival time of the process; 0 if the total time is censored and 1 otherwise.}
  \item{\dots}{Other options.}
}

\value{An object of class "survCS" and of class "surv". "survCS" objected are implemented as a single element list. }

\author{Luis Meira-Machado and Marta Sestelo}

\examples{
data(colonCS)
obj <- with(colonCS, survCS(time1, event1, Stime, event, age = age))
obj1 <- with(colonCS, survCS(time1, event1, Stime, event, age = age, nodes = nodes))
}

