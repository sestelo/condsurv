#' Kaplan-Meier weights.
#' 
#' This function returns a vector with the Kaplan-Meier weights.
#' 
#' 
#' @param time Survival time of the process.
#' @param status Censoring indicator of the survival time of the process; 0 if
#' the survival time is censored and 1 otherwise.
#' @return Vector with Kaplan-Meier weights.
#' @author Luis Meira-Machado and Marta Sestelo
#' @seealso \code{\link{PKMW}}
#' @references E. Kaplan and P. Meier. Nonparametric estimation from incomplete
#' observations. Journal of the American Statistical Association, 53:457-481,
#' 1958.
#' @examples
#' 
#' obj <- with(colonCS, survCS(time1, event1, Stime, event))
#' kmw <- KMW(time = obj$Stime, status = obj$event)
#' 
#' require(survival)
#' colon.surv <- survfit(Surv(Stime, event) ~ 1, obj)
#' times <- summary(colon.surv)$time
#' surv <- summary(colon.surv)$surv
#' nevent <- summary(colon.surv)$n.event
#' p <- match(obj$Stime, times)
#' kmw2 <- -diff(c(1, surv))/nevent
#' kmw2 <- kmw2[p]*obj$event
#' kmw2[is.na(kmw2)] <- 0
#' all.equal(kmw, kmw2)
#' 
#' 
KMW <-
  function(time, status){
    t1 <- max(time)
    len <- length(time)
    res <- .C( "WeightsKaplanMeierSort", time = as.double(time),
               status = as.integer(status), as.integer(len),
               as.double(t1), weights = double(len) )
    return(res$weights)
  }
