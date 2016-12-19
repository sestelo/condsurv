\name{survCS}
\alias{survCS}
\title{Create a survCS object.}
\description{Creates a "survCS" object, usually used as a response variable in a model formula.}
\usage{
survCS(time1, event1, Stime, event, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time1}{First time or censoring time.}
  \item{event1}{Indicator of the first time; 0 if the first time is censored and 1 otherwise.}
  \item{Stime}{The total time of the process.}
  \item{event}{Censoring indicator of the survival time of the process;
  0 if the total time is censored and 1 otherwise.}
  \item{\dots}{Other options.}
}

\details{ Arguments in this function must be introduced in the following
order: time1, event1, time2, ..., Stime and event, where time1,
time2, ..., Stime are ordered event times and event1, event2,..., event
their corresponding indicator statuses.}

\value{An object of class "survCS" and of class "surv". "survCS" objects
are implemented as a single dataframe. }

\author{Luis Meira-Machado and Marta Sestelo}

\examples{

with(colonCS, survCS(time1, event1, Stime, event))

}

