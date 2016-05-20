\name{survIPCW}
\alias{survIPCW}

\title{Inverse probability of censoring weighted conditional survival probabilities.
}
\description{Provides estimates for the conditional survival probabilities based on inverse probability of censoring weighted estimators, conditionally on current or past covariate measures.}
\usage{
survIPCW(object, x, y, z.name, z.value, bw = "dpik", window = "gaussian",
method.weights = "NW", conf = FALSE, n.boot = 200, conf.level = 0.95,
lower.tail = FALSE, cluster = TRUE, ncores = NULL)
}

\arguments{
  \item{object}{An object of class "survCS".}
  \item{x}{The first time for obtaining estimates for the conditional survival probabilities.}
  \item{y}{The total time for obtaining estimates for the conditional survival probabilities.}
  \item{z.name}{The name of the covariate.}
  \item{z.value}{The value of the covariate ("z.name") at which the conditional survival probabilities are computed.}
  \item{bw}{A single numeric value to compute a kernel density bandwidth.}
  \item{window}{A character string specifying the desired kernel. See details below for possible options. Defaults to "gaussian" where the gaussian density kernel will be used.}
  \item{method.weights}{A character string specifying the desired weights method. Possible options are "NW" for the Nadaraya-Watson weights and "LL" for local linear weights. Defaults to "NW".}
  \item{conf}{Provides pointwise confidence bands. Defaults to FALSE.}
  \item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
  \item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
  \item{lower.tail}{logical; if FALSE (default), probabilities are P(T > y|T1 > x) otherwise, P(T > y|T1 <= x).}
  \item{cluster}{A logical value. If  \code{TRUE} (default), the
bootstrap procedure for the confidence intervals is  parallelized. Note that there are cases
(e.g., a low number of bootstrap repetitions) that R will gain in
performance through serial computation. R takes time to distribute tasks
across the processors also it will need time for binding them all together
later on. Therefore, if the time for distributing and gathering pieces
together is greater than the time need for single-thread computing, it does
not worth parallelize.}

\item{ncores}{An integer value specifying the number of cores to be used
in the parallelized procedure. If \code{NULL} (default), the number of cores
to be used is equal to the number of cores of the machine - 1.}

}


\details{Possible options for argument window are "gaussian", "epanechnikov", "tricube", "boxcar", "triangular", "quartic" or "cosine".}


\value{An object of class "IPCW". "IPCW" objects are implemented as a list with elements:
  \item{est}{data.frame with estimates of the conditional probabilities.}
  \item{estimate}{Estimates of the conditional survival probability.}
  \item{LCI}{The lower conditional survival probabilities of the interval.}
  \item{UCI}{The upper conditional survival probabilities of the interval.}
  \item{conf.level}{Level of confidence.}
  \item{z.name}{The name of the covariate.}
  \item{z.value}{The value of the covariate.}
\item{bw}{A character string indicating a function to compute
a kernel density.}
\item{window}{A character string specifying the desired kernel.}
\item{method.weights}{A character string specifying
the desired weights method.}
  \item{x}{The first time for obtaining the estimates of the conditional survival probabilities.}
  \item{y}{The total time for obtaining the estimates of the conditional survival probabilities.}
  \item{conf}{logical; if FALSE (default) the pointwise confidence bands are not given.}
   \item{conf}{logical; if FALSE (default) the pointwise confidence bands are not given.}
    \item{lbd}{bandwidth.}
}


\references{
L. Meira-Machado, M. Sestelo, and A. Goncalves. Nonparametric estimation of the survival function
for ordered multivariate failure time data: a comparative study. Biometrical Journal, 2016.
}

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{survKMW}}, \code{\link{survLDM}} and \code{\link{survPLDM}}}

\examples{

data(colonCS)
obj2 <- with(colonCS, survCS(time1, event1, Stime, event, age = age))
survLDM(obj2, x = 365, y = 730)
survIPCW(obj2, x = 365, y = 730, z.name = "age", z.value = 48)
survIPCW(obj2, x = 365, y = 730, z.name = "age",
z.value = c(30, 35, 40, 45, 50, 55, 60))
survIPCW(obj2, x = 365, y = c(730, 1095, 1460),
z.name = "age", z.value = c(48))
}

