survCOND <-
  function(object, x, y, lower.tail = FALSE, method = "LDM",
           presmooth = FALSE, conf = FALSE, n.boot = 1000,
           conf.level = 0.95, cluster = FALSE, ncores = NULL, na.rm = TRUE)
  {

      if (!(method %in% c("LDM", "KMW"))){
        stop("Possible methods are 'LDM' and 'KMW'." )
        }

    if (method == "KMW" & presmooth == TRUE)  {
      warning("Argument 'presmooth' is not used by the 'KMW' method.")
    }

    if (method == "LDM" & presmooth == FALSE)  {
      res <- survLDM(object = object, x = x, y = y, conf = conf,
                     n.boot = n.boot, conf.level = conf.level,
                     lower.tail = lower.tail, cluster = cluster,
                     ncores = ncores)

      class(res) <- c("LDM", "surv")
      return(invisible(res))
    }

    if (method == "LDM" & presmooth == TRUE) {
      res <- suppressWarnings(survPLDM(object = object, x = x, y = y,
                                       conf = conf, n.boot = n.boot,
                                       conf.level = conf.level,
                                       lower.tail = lower.tail,
                                       cluster = cluster, ncores = ncores))

      class(res) <- c("PLDM", "surv")
      return(invisible(res))
    }

    if (method == "KMW")  {
      res <- survKMW(object = object, x = x, y = y, conf = conf,
                     n.boot = n.boot, conf.level = conf.level,
                     lower.tail = lower.tail, cluster = cluster,
                     ncores = ncores, na.rm = na.rm)

      class(res) <- c("KMW", "surv")
      return(invisible(res))
    }
  }

