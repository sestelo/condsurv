\name{survCOND}
\alias{survCOND}

\title{Conditional survival probabilities based on the Kaplan-Meier weights and Landmark approaches.}
\description{Provides estimates for the conditional survival probabilities
  based on Kaplan-Meier weighted estimators and the Landmark approaches.}
% \usage{
%   survCOND(object, x, y, lower.tail = FALSE, method = "LDM", presmooth = FALSE, conf = FALSE, n.boot = 1000, conf.level = 0.95,
%           cluster = FALSE, ncores = NULL, na.rm = TRUE)
% }
%
% \arguments{
%   \item{object}{An object of class "survCS".}
%   \item{x}{Time or vector of times for the condional event(s).}
%   \item{y}{The total time for obtaining estimates for the conditional survival probabilities.}
%   \item{lower.tail}{vector of logical values with the same size as 'x'. If 'x' has dimension
%     one and if lower.tail = FALSE (default), probabilities are P(T > y|T1 > x) otherwise, P(T > y|T1 <= x).
%     If the conditional event is 2-dimensional, then, for example, given x = c(x1, x2) and lower.tail = c(TRUE, FALSE)
%     must be used to obtain probabilities P(T > y|T1 <= x1, T2 > x2). Multi-dimensional conditional events are introduced similarly.}
%   \item{method}{The method used to compute the conditional survival function.
%     Possible options are "LDM" and "KMW". Defaults to "LDM".
%   \item{presmooth}{A logical value. If  \code{TRUE}, the presmoothed landmark estimator of the
%     conditional survival function is computed.
%   \item{conf}{Provides pointwise confidence bands. Defaults to FALSE.}
%   \item{n.boot}{The number of bootstrap samples. Defaults to 1000 samples.}
%   \item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}
%   \item{cluster}{A logical value. If  \code{TRUE} (default), the
%     bootstrap procedure  for the confidence intervals is  parallelized. Note that there are cases
%     (e.g., a low number of bootstrap repetitions) that R will gain in
%     performance through serial computation. R takes time to distribute tasks
%     across the processors also it will need time for binding them all together
%     later on. Therefore, if the time for distributing and gathering pieces
%     together is greater than the time need for single-thread computing, it does
%     not worth parallelize.}
%
%   \item{ncores}{An integer value specifying the number of cores to be used
%     in the parallelized procedure. If \code{NULL} (default), the number of cores
%     to be used is equal to the number of cores of the machine - 1.}
%   \item{na.rm}{A logical value indicating whether NA values should be stripped before the computation proceeds.}
%
% }
%
% \value{An object of class "surv" and one of the following three classes: "KMW", "LMD" and "PLDM". Objects are implemented as a list with elements:
%     \item{est}{data.frame with estimates of the conditional probabilities.}
%   \item{estimate}{Estimates of the conditional survival probability.}
%   \item{LCI}{The lower conditional survival probabilities of the interval.}
%   \item{UCI}{The upper conditional survival probabilities of the interval.}
%   \item{conf.level}{Level of confidence.}
%   \item{y}{The total time for obtaining the estimates of the conditional survival probabilities.}
%   \item{x}{The first time for obtaining the estimates of the conditional survival probabilities.}
%   \item{conf}{logical; if FALSE (default) the pointwise confidence bands are not given.}
% }
%
% \references{
%   L. Meira-Machado, M. Sestelo, and A. Goncalves. Nonparametric estimation of the survival function
%   for ordered multivariate failure time data: a comparative study. Biometrical Journal, 2016.
% }
%
% \author{Luis Meira-Machado and Marta Sestelo}
%
%
% \seealso{\code{\link{survIPCW}}}
%
% \examples{
%   data(colonCS)
%   obj <- with(colonCS, survCS(time1, event1, Stime, event))
%   %P(T>y|T1>x)
%   survCOND(obj, x = 365, y = 730, method = "KMW")
%
%   survCOND(obj, x = 365, y = c(730, 1095, 1460), method = "KMW")
%
%   %P(T>y|T1<=x)
%   survCOND(obj, x = 365, y = 730, method = "KMW", lower.tail = TRUE)
%
%   survCOND(obj, x = 365, y = 730, method = "LMD, lower.tail = TRUE)
%
%   survCOND(obj, x = 365, y = 730, method = "LMD, presmooth = TRUE,
%       lower.tail = TRUE, conf = TRUE, n.boot = 100, conf.level = 0.95, cluster = FALSE)
%
%   res <- survCOND(obj, x = 365, y = c(730, 1095, 1460), method = "LMD,
%       lower.tail = TRUE, conf = TRUE, n.boot = 100, conf.level = 0.95, cluster = FALSE)
%   res$est
% }
