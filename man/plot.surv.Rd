\name{plot.surv}
\alias{plot.surv}

\title{Plot for an object of class "surv".}
\description{It draws the estimated conditional survival probabilities. }
\usage{
\method{plot}{surv}(x = object, y = NULL, type = NULL, ...)
}

\arguments{
  \item{x}{An object of class "surv".}
  \item{y}{NULL}
  \item{type}{The type of plot should be drawn. See details below for possible options. Defaults to "l" for the draw be a line.}
  \item{\dots}{Other options.}
}

\value{No value is returned. }

\author{Luis Meira-Machado and Marta Sestelo}


\examples{
data(colonCS)
obj <- with(colonCS, survCS(time1, event1, Stime, event, age = age))
y1 <- obj$data$Stime
y1 <- y1[y1 > 365]
res1 <- survLDM(obj, x = 365, y = y1, conf = TRUE, n.boot = 100,
conf.level = 0.95, cluster = FALSE)
plot(res1, type = "l", col = c(1, 2, 2), lty = c(1, 2, 2),
xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))
}

