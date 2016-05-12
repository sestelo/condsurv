\name{survLDM}
\alias{survLDM}

\title{Conditional survival probabilities based on the Landmark approach.}
\description{Provides estimates for the conditional probabilities based on the landmark approach. 
Given the time point x, to estimate $P(T > y | T1 > x)$ the analysis can be restricted to the individuals with an observed
first event time greater than x.}

\usage{
survLDM(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95, lower.tail = FALSE)
}

\arguments{
  \item{object}{An object of class "survCS".}
  \item{x}{The first time for obtaining estimates for the conditional survival probabilities.}
  \item{y}{The total time for obtaining estimates for the conditional survival probabilities.}
  \item{conf}{Provides pointwise confidence bands. Defaults to FALSE.}
  \item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
  \item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
  \item{lower.tail}{logical; if FALSE (default), probabilities are P(T > y|T1 > x) otherwise, P(T <= y|T1 <= x).}
}


\value{An object of class "LDM" and of the class "surv". "LDM" objects are implemented as a list with elements: 
  \item{est}{data.frame with estimates of the conditional probabilities.}
  \item{estimate}{Estimates of the conditional survival probability.}
  \item{LCI}{The lower conditional survival probabilities of the interval.}
  \item{UCI}{The upper conditional survival probabilities of the interval.}
  \item{conf.level}{Level of confidence.}
  \item{y}{The total time for obtaining the estimates of the conditional survival probabilities.}
  \item{x}{The first time for obtaining the estimates of the conditional survival probabilities.}
  \item{conf}{logical; if FALSE (default) the pointwise confidence bands are not given.}
}

\references{
L. Meira-Machado, M. Sestelo, and A. Goncalves. Nonparametric estimation of the survival function
for ordered multivariate failure time data: a comparative study. Biometrical Journal, 2016.

H. Van Houwelingen. Dynamic prediction by landmarking in event history analysis. Scandinavian
Journal of Statistics, 34:70-85, 2007.
}

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{survKMW}}, \code{\link{survPLDM}} and \code{\link{survIPCW}}}

\examples{
data(colonCS)
obj <- with(colonCS, survCS(time1, event1, Stime, event))

#P(T > y|T > x)
survLDM(obj, x = 365, y = 730)

#P(T > y|T <= x)
survLDM(obj, x = 365, y = 730, lower.tail = TRUE)
survLDM(obj, x = 365, y = c(730, 1095, 1460))  
survLDM(obj, x = 365, y = 730, conf = TRUE, n.boot = 500, conf.level = 0.95)
survLDM(obj, x = 365, y = c(730, 1095, 1460), conf = TRUE, n.boot = 500, conf.level = 0.95)

res <- survLDM(obj, x = 365, y = c(730, 1095, 1460), conf = TRUE, n.boot = 500, conf.level = 0.95)
res$est
}
