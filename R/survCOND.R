survCOND <-
  function(object, x, y, lower.tail = FALSE, method = "LDM", conf = FALSE,
           n.boot = 1000, conf.level = 0.95, cluster = FALSE,
           ncores = NULL, na.rm=T)
{

   if (method == "LDM")  {
	   res <- survLDM(object = object, x = x, y = y, conf = conf,
	   n.boot = n.boot, conf.level = conf.level, lower.tail = lower.tail,
           cluster = cluster, ncores = ncores)

	   class(res) <- c("LDM", "surv")
           return(invisible(res))
				  }

   if (method == "PLDM") {
	   res <- suppressWarnings(survPLDM(object = object, x = x, y = y, conf = conf,
	   n.boot = n.boot, conf.level = conf.level, lower.tail = lower.tail,
           cluster = cluster, ncores = ncores))

	   class(res) <- c("PLDM", "surv")
           return(invisible(res))
				  }

   if (method == "KMW")  {
	   res <- survKMW(object = object, x = x, y = y, conf = conf,
	   n.boot = n.boot, conf.level = conf.level, lower.tail = lower.tail,
           cluster = cluster, ncores = ncores, na.rm=na.rm)

	   class(res) <- c("KMW", "surv")
           return(invisible(res))
				  }
}


