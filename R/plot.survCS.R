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



      if (is.null(ylim)) {
        if (object$Nlevels == 1) {
          ymin <- ifelse(ci == TRUE, min(ob[, 3]), min(ob[, 2]))
        }else{
          ymin <- ifelse(ci == TRUE,
                         min(sapply(ob, function(x) min(x[, 3]), simplify = TRUE)),
                         min(sapply(ob, function(x) min(x[, 2]), simplify = TRUE)))
        }
        ylim <- c(ymin, 1)
      }


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






