#' Presmoothed Kaplan-Meier weights.
#' 
#' This function returns a vector with the presmoothed Kaplan-Meier weights.
#' 
#' 
#' @param time Survival time of the process.
#' @param status Censoring indicator of the survival time of the process; 0 if
#' the survival time is censored and 1 otherwise.
#' @return Vector with presmoothed Kaplan-Meier weights.
#' @author Luis Meira-Machado and Marta Sestelo
#' @seealso \code{\link{KMW}}
#' @references R. Cao, I. Lopez-de Ullibarri, P. Janssen, and N. Veraverbeke.
#' Presmoothed kaplan-meier and nelsonaalen estimators. Journal of
#' Nonparametric Statistics, 17:31-56, 2005.
#' 
#' G. Dikta. On semiparametric random censorship models. Journal of Statistical
#' Planning and Inference, 66:253-279, 1998.
#' 
#' E. Kaplan and P. Meier. Nonparametric estimation from incomplete
#' observations. Journal of the American Statistical Association, 53:457-481,
#' 1958.
#' @examples
#' 
#' obj <- with(colonCS, survCS(time1, event1, Stime, event))
#' PKMW(time = obj$Stime, status = obj$event)
#' 
PKMW <-
  function(time, status){
    t1 <- max(time)
    len <- length(time)
    status1 <- fitted(glm(status ~ time, family = binomial))
    res <- .C( "WeightsKaplanMeierSortEx", time=as.double(time),
               status <- as.double(status1), as.integer(len), as.double(t1),
               weights = double(len) )
    return(res$weights)
  }
