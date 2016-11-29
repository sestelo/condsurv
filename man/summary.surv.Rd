\name{summary.surv}
\alias{summary.surv}
\alias{print.surv}
\title{Summarizing fits of \code{surv} class}
\description{Returns a a data.frame or list containing the estimates of
the conditional survival, its confidence limits and other information. }
\usage{
  \method{summary}{surv}(object, times = NULL, ...)
}

  \arguments{
\item{object}{A fitted \code{surv} object as produced by \code{survCOND()}.}
\item{times}{Vector of times; the returned data frame will contain 1 row
for each time. Missing values are not allowed.}
\item{...}{For future methods.}
}

\value{A data frame or a list containing the following components:
 \item{y}{The total time for obtaining the estimates of the conditional survival probabilities.}
  \item{est}{Estimates of the conditional survival probability.}
  \item{lower 95\% CI}{The lower conditional survival probabilities of the interval.}
  \item{upper 95\% CI}{The upper conditional survival probabilities of the interval.}
}

    \examples{
    fit <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
    data = colonCS, method = "LDM", conf = TRUE, n.boot = 50, cluster = FALSE)
   summary(fit)
    summary(fit, times = c(400, 1000, 2900))

    }

\references{
L. Meira-Machado, M. Sestelo, and A. Goncalves (2016). Nonparametric estimation of the survival function
for ordered multivariate failure time data: a comparative study. Biometrical Journal, 58(3), 623--634.
}

   \author{Luis Meira-Machado and Marta Sestelo}



