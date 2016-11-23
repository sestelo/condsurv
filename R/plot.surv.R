plot.surv <-
  function(x = object, y = NULL, conf = NULL,  type = NULL, conftype = NULL,
           col = 1, confcol = 2, lty = 1, conflty = 2, ...) {

    object <- x

    if (class(object)[1] != "KMW" & class(object)[1] != "LDM" &
        class(object)[1] != "PLDM" &
        class(object)[1] != "IPCW"){

      stop("The argumment 'Object' must be of one of the following classes
         'KMW', 'LDM', 'PLDM', 'IPCW'")
    }

    ob <- object$est

    if (is.null(type)) type <- "l"
    if (is.null(conftype)) conftype <- "l"


    if (is.null(conf)){
      ci <- object$conf
    }else{
      if (conf == TRUE & object$conf == FALSE){
        stop("The surv object does not contain confidence intervals")
      }
      if (conf == FALSE) ci <- FALSE
    }


      plot(ob[, 1], ob[, 2], type = type, col = col,...)
      if (ci == TRUE){
        lines(x = ob[, 1], y = ob[, 3], type = conftype, lty = conflty,
              col = confcol, ...)
        lines(x = ob[, 1], y = ob[, 4], type = conftype, lty = conflty,
              col = confcol, ...)
      }


    #
    #
    # if (is.null(conf) & object$conf == FALSE) { plot(ob[, 1], ob[, 2], type = type, ...) }
    #
    # if (is.null(conf) & object$conf == TRUE) {
    #   matplot(x = ob[, 1], y = ob[, 2:4], type = type, ...)
    #   #oask <- devAskNewPage(TRUE)
    #   #on.exit(devAskNewPage(oask))
    # }

  }
