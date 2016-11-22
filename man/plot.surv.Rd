\name{plot.surv}
\alias{plot.surv}

\title{Plot for an object of class "surv".}
\description{It draws the estimated conditional survival probabilities. }
\usage{
\method{plot}{surv}(x = object, y = NULL, conf = NULL,  type = NULL, conftype = NULL,
           col = 1, confcol = 2, lty = 1, conflty = 2, ...)
}

\arguments{
  \item{x}{An object of class "surv".}
  \item{y}{\code{NULL}}
  \item{conf}{Draw the confidence intervals into the plot. By default it is
  \code{NULL} (they are drawn if the "surv" object contains them.}
  \item{type}{The type of plot should be drawn. See details \code{\link{par}}
  for possible options. Defaults to "l" for the draw be a line.}
  \item{conftype}{The type of plot should be drawn for confidence intervals. See details \code{\link{par}}
  for possible options. Defaults to "l" for the draw be a line.}
  \item{col}{A specification for the default plotting color.}
  \item{confcol}{A specification for the default confidence intervals
  plotting color}
  \item{lty}{The line type. Line types can either be specified as an integer
  (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash,
  5 = longdash, 6 = twodash).  See details in \code{\link{par}}.}
  \item{conflty}{The line type for confidence intervals. Line types can either
  be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
  3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).}
  \item{\dots}{Other options.}
}

\value{No value is returned. }

\author{Luis Meira-Machado and Marta Sestelo}


\examples{
data(colonCS)
obj <- with(colonCS, survCS(time1, event1, Stime, event, age = age))

res1 <- survCOND(obj, x = 365, y = obj$data$Stime, method = "LDM",
conf = TRUE, n.boot = 100, conf.level = 0.95, cluster = FALSE)

plot(res1, type = "l", xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))

plot(res1, type = "l", conf = FALSE, xlab = "Time (days)",
ylab = "S(y|365)", ylim = c(0.5, 1))
}

