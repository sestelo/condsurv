#' Kaplan-Meier product-limit estimate of survival.
#' 
#' This function provides survival estimates using the product-limit
#' Kaplan-Meier estimator.
#' 
#' 
#' @param time Survival time of the process.
#' @param status Censoring indicator of the survival time of the process; 0 if
#' the survival time is censored and 1 otherwise.
#' @param t The time for obtaining survival estimates.
#' @author Luis Meira-Machado and Marta Sestelo
#' @seealso \code{\link{PKM}}
#' @references E. Kaplan and P. Meier. Nonparametric estimation from incomplete
#' observations. Journal of the American Statistical Association, 53:457-481,
#' 1958.
#' @examples
#' 
#' require(survival)
#' obj <- with(colonCS, survCS(time1, event1, Stime, event))
#' KM(time = obj$Stime, status = obj$event, t = 1095)
#' 
#' fit <- survfit(Surv(obj$Stime, obj$event) ~ 1, data = obj)
#' summary(fit, time = 1095)$surv
#' 
#' 
KM <-
  function(time, status, t)
  {
    if (missing(time)) stop("Argument 'time' is missing with no default")
    if (missing(status)) stop("Argument 'status' is missing with no default")
    if (missing(t)) stop("Argument 't' is missing with no default")
    return( .C( "KaplanMeierValueSort", as.double(time), as.integer(status),
                as.integer( length(time) ), as.double(t), as.double(1) )[[5]] )
  }
