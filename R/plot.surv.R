plot.surv <- function(x = object, y = NULL, conf = NULL, type = NULL,
  conftype = NULL, col = 1:6, confcol = 1:6, lty = 1, conflty = 2,
  xlab = "Time", ylab = "Probability", ...) {

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


  if (object$Nlevels == 1) {
    plot(ob[, 1], ob[, 2], type = type, col = col, xlab = xlab,
      ylab = ylab, lty = lty, ...)
    if (ci == TRUE) {
      lines(x = ob[, 1], y = ob[, 3], type = conftype,
        lty = conflty, col = confcol, ...)
      lines(x = ob[, 1], y = ob[, 4], type = conftype,
        lty = conflty, col = confcol, ...)
    }
  } else {

    plot(ob[[1]][, 1], ob[[1]][, 2], type = "n", xlab = xlab,
      ylab = ylab, ...)

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






