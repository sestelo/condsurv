#' Plot for an object of class "survCS".
#' 
#' It draws the estimated conditional survival probabilities.
#' 
#' 
#' @param x An object of class "survCS".
#' @param y \code{NULL}
#' @param conf Draw the confidence intervals into the plot. By default it is
#' \code{NULL}, they are drawn if the "surv" object contains them.
#' @param type The type of plot that should be drawn. See details
#' \code{\link{par}} for possible options. Defaults to "s" for the draw be
#' stair steps.
#' @param conftype The type of plot that should be drawn for confidence
#' intervals.  See details \code{\link{par}} for possible options. Defaults to
#' "s" for the draw be stair steps.
#' @param col Vector of colors. Colors are used cyclically.
#' @param confcol Vector of colors for the confidence intervals.  Colors are
#' used cyclically.
#' @param lty The line type. Line types can either be specified as an integer
#' (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash, 5 =
#' longdash, 6 = twodash).  See details in \code{\link{par}}.
#' @param conflty The line type for confidence intervals. Line types can either
#' be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed, 3 =
#' dotted, 4 = dotdash, 5 = longdash, 6 = twodash).
#' @param xlab A title for the \code{x} axis: see \code{\link{title}}.
#' @param ylab A title for the \code{y} axis: see \code{\link{title}}.
#' @param ylim The \code{y} limits of the plot.
#' @param xlim The \code{x} limits of the plot.
#' @param \dots Other options.
#' @return No value is returned.
#' @author Luis Meira-Machado and Marta Sestelo
#' @examples
#' 
#' fit1 <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
#'    data = colonCS, method = "LDM", conf = TRUE)
#' 
#' plot(fit1, xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))
#' 
#' fit4 <- survCOND(survCS(time1, event1, Stime, event) ~ rx,
#'    x = 365, data = colonCS, method = "LDM")
#' 
#' plot(fit4, xlab = "Time (days)", ylab = "S(y|365)", ylim = c(0.5, 1))
#' 
#' 
#' @export plot.survCS
plot.survCS <- function(x = object, y = NULL, conf = NULL, type = NULL,
                        conftype = NULL, col = 1:6, confcol = 1:6, lty = 1, conflty = 2,
                        xlab = "Time", ylab = "Survival", ylim = NULL, xlim = NULL, ...) {

  if (inherits(x, "survCS")) {

    if (class(x)[1] == "data.frame") {
      plot(x)

    }else{

      object <- x

      if (object$Nlevels != length(col))
        col <- rep(col, times = object$Nlevels)
      if (object$Nlevels != length(confcol))
        confcol <- rep(confcol, times = object$Nlevels)

      if (class(object)[1] != "KMW" & class(object)[1] != "LDM" &
          class(object)[1] != "PLDM" & class(object)[1] != "IPCW") {
        stop("The argumment 'Object' must be of one of the following classes
     'KMW', 'LDM', 'PLDM', 'IPCW'")
      }

      ob <- object$est

      if (is.null(type))
        type <- "s"
      if (is.null(conftype))
        conftype <- "s"


      if (is.null(conf)) {
        ci <- object$conf
      } else {
        if (conf == TRUE & object$conf == FALSE) {
          stop("The surv object does not contain confidence intervals")
        }
        if (conf == TRUE & object$conf == TRUE)
          ci <- TRUE
        if (conf == FALSE)
          ci <- FALSE

      }


#
#       if (is.null(ylim)) {
#         if (object$Nlevels == 1) {
#           ymin <- ifelse(ci == TRUE, min(ob[, 3]), min(ob[, 2]))
#         }else{
#           ymin <- ifelse(ci == TRUE,
#                          min(sapply(ob, function(x) min(x[, 3]), simplify = TRUE)),
#                          min(sapply(ob, function(x) min(x[, 2]), simplify = TRUE)))
#         }
#         ylim <- c(ymin, 1)
#       }


      if (is.null(ylim)) ylim <- c(0, 1)


      if (is.null(xlim) & object$Nlevels > 1) {
        xlim <- c(min(sapply(ob, function(x) min(x[, 1]), simplify = TRUE)),
                  max(sapply(ob, function(x) max(x[, 1]), simplify = TRUE)))
      }


      #if (is.null(ylim)) ylim <- c(0,1)

      if (object$Nlevels == 1) {

        plot(ob[, 1], ob[, 2], type = type, col = col, xlab = xlab,
             ylab = ylab, lty = lty, ylim = ylim, xlim = xlim, ...)
        if (ci == TRUE) {
          lines(x = ob[, 1], y = ob[, 3], type = conftype,
                lty = conflty, col = confcol, ...)
          lines(x = ob[, 1], y = ob[, 4], type = conftype,
                lty = conflty, col = confcol, ...)
        }
      } else {

        plot(ob[[1]][, 1], ob[[1]][, 2], type = "n", xlab = xlab,
             ylab = ylab, ylim = ylim, xlim = xlim, ...)

        for (i in 1:object$Nlevels) {
          lines(ob[[i]][, 1], ob[[i]][, 2], type = type, col = col[i],
                lty = lty, ...)
          if (ci == TRUE) {
            lines(x = ob[[i]][, 1], y = ob[[i]][, 3], type = conftype,
                  lty = conflty, col = confcol[i], ...)
            lines(x = ob[[i]][, 1], y = ob[[i]][, 4], type = conftype,
                  lty = conflty, col = confcol[i], ...)

          }
        }
        legend("topright", object$levels, col = col, lty = lty)
      }
    }
  }else{
    stop("Argument x must be either survCS object.")
  }
}






