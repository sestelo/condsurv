#' Summarizing fits of \code{survCS} class
#'
#' Returns a a data.frame or list containing the estimates of the conditional
#' survival, its confidence limits and other information.
#'
#' @aliases summary.survCS print.survCS
#' @param object A fitted \code{survCS} object as produced by
#' \code{survCOND()}.
#' @param times Vector of times; the returned data frame will contain 1 row for
#' each time. Missing values are not allowed.
#' @param ... For future methods.
#'
#' @return A data frame or a list containing the following components:
#' \item{y}{The total time for obtaining the estimates of the conditional
#' survival probabilities.}
#' \item{est}{Estimates of the conditional survival probability.}
#' \item{lower 95\% CI}{The lower conditional survival probabilities of the interval.}
#' \item{upper 95\% CI}{The upper conditional survival probabilities of the interval.}
#'
#' @author Luis Meira-Machado and Marta Sestelo
#'
#' @references L. Meira-Machado, M. Sestelo, and A. Goncalves (2016).
#' Nonparametric estimation of the survival function for ordered multivariate
#' failure time data: a comparative study. Biometrical Journal, 58(3),
#' 623--634.
#'
#' @examples
#' fit <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
#' data = colonCS, method = "LDM", conf = TRUE, n.boot = 50, cluster = FALSE)
#' summary(fit)
#' summary(fit, times = c(400, 1000, 2900))
#'
#'
#' @export summary.survCS
summary.survCS <- function(object, times = NULL, ...){


  if (inherits(object, "survCS")) {

    if (class(object)[1] == "data.frame"){
      summary(object)
    }else{

      if (is.null(times)) {
        cat("\n")
        cat(object$callp, "\n")
        cat("\n")
        if (object$Nlevels > 1) {
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            print(object$est[[i]], row.names = FALSE)
            cat("\n")
          }
        }else{
          print(object$est, row.names = FALSE)
        }
        res <- object$est
      }else{
        if (object$Nlevels > 1) {

          # to control the times argument
          # -----------------------------
          pp <- list()
          for (i in 1:object$Nlevels) {
            pp[[i]] <- sapply(times, function(x)ifelse(x >= min(object$est[[i]]$y) & x <= max(object$est[[i]]$y), 1, NA))
          }
          if (all(is.na(unlist(pp)))) {
            stop(paste("At least one element of the 'times' vector has to be in the range of 'y' "))
          }

          if (any(is.na(unlist(pp)))) {
            warning(paste("'times' must be in the range of 'y' (for each level)" ))
          }
          #--------------------

          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            ii <- sapply(times, function(x)ifelse( x >= min(object$est[[i]]$y) & x <= max(object$est[[i]]$y),
                                                   which.max(object$est[[i]]$y[object$est[[i]]$y <= x]), NA))

            if (all(is.na(ii))) {
              aux <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[i]])[2] - 1))
            }else{
              aux <- data.frame(times, object$est[[i]][ii,-1 ])
            }
            names(aux) <- names(object$est[[i]])
            res[[paste(object$levels[i])]] <- aux
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            print(res[[i]], row.names = FALSE)
            cat("\n")
          }
        }else{
          ii <- sapply(times, function(x)ifelse( x >= min(object$est$y) & x <= max(object$est$y),
                                                 which.max(object$est$y[object$est$y <= x]), NA))
          if (all(is.na(ii))) {
            stop(paste("At least one element of the 'times' vector has to be between",min(object$est$y), "and", max(object$est$y)))
          }

          if (any(is.na(ii))) {
            warning(paste("'times' must be between",min(object$est$y), "and", max(object$est$y)))
          }
          res <- data.frame(times, object$est[ii, -1])
          names(res) <- names(object$est)
          print(res, row.names = FALSE)
        }
      }

      class(res) <- "summary.surv"
      return(invisible(res))

    }
  }else{
    stop("Argument x must be either survCS object.")
  }
}


