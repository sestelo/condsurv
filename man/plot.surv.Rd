\name{plot.surv}
\alias{plot.surv}
\title{Plot for an object of class "surv".}
\description{It draws the estimated conditional survival probabilities. }
\usage{
\method{plot}{surv}(x = object, y = NULL, conf = NULL,  type = NULL, conftype = NULL,
           col = 1:6, confcol = 1:6, lty = 1, conflty = 2,
           xlab = "Time", ylab = "Probability", ...)
}

\arguments{
  \item{x}{An object of class "surv".}
  \item{y}{\code{NULL}}
  \item{conf}{Draw the confidence intervals into the plot. By default it is
  \code{NULL} (they are drawn if the "surv" object contains them.}
  \item{type}{The type of plot should be drawn. See details \code{\link{par}}
  for possible options. Defaults to "l" for the draw be a line.}
  \item{conftype}{The type of plot should be drawn for confidence intervals.
  See details \code{\link{par}}
  for possible options. Defaults to "l" for the draw be a line.}
  \item{col}{Vector of colors. Colors are used cyclically.}
  \item{confcol}{Vector of colors for the confidence intervals.
  Colors are used cyclically.}
  \item{lty}{The line type. Line types can either be specified as an integer
  (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash,
  5 = longdash, 6 = twodash).  See details in \code{\link{par}}.}
  \item{conflty}{The line type for confidence intervals. Line types can either
  be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
  3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).}
  \item{xlab}{A title for the \code{x} axis: see \code{\link{title}}.}
  \item{ylab}{A title for the \code{y} axis: see \code{\link{title}}.}
  \item{\dots}{Other options.}
}

\value{No value is returned. }

\author{Luis Meira-Machado and Marta Sestelo}


\examples{
fit1 <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
   data = colonCS, method = "LDM", conf = TRUE)

plot(fit1, xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))

fit4 <- survCOND(survCS(time1, event1, Stime, event) ~ rx,
   x = 365, data = colonCS, method = "LDM")

plot(fit4, xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))

}

