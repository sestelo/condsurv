survLDM <-
  function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95,
           lower.tail = FALSE, cluster = FALSE, ncores = NULL)
  {
    lenc <- dim(object[[1]])[2]
    ntimes <- lenc%/%2

    X <- data.frame(object[[1]][,2*(1:ntimes)-1])
    p1 <- whichCS(X, x=x, lower.tail=lower.tail)
    if (length(p1) == 0) stop("Insufficient data.")
    G0 <- KMW(object[[1]]$Stime[p1], object[[1]]$event[p1])
    t_2 <- object[[1]]$Stime[p1]
    res <- rep(0, length(y))
    for (k in 1: length(y)) {
      p2 <- which(t_2 <= y[k])
      res[k] <- 1 - sum(G0[p2])
      if (res[k] < 0) res[k] <- 0
    }


    ii <- duplicated(res)
    y <- y[!ii]

    res.li <- rep(0, length(y))
    res.ls <- rep(0, length(y))
    resu <- data.frame(cbind(y, res[!ii]))
    names(resu) <- c("y", "estimate")


    if (conf==TRUE) {
      simplebootsurvLDM <- function(object, y, x, lower.tail){
        j <- 1
        res.ci <- matrix(0, nrow = length(y), ncol = j)
        n <- dim(object[[1]])[1]
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        X <- data.frame(ndata[,2*(1:ntimes)-1])
        p1 <- whichCS(X, x=x, lower.tail=lower.tail)
        if (length(p1) == 0) stop("Insufficient data.")
        G0 <- KMW(ndata$Stime[p1], ndata$event[p1])
        t_2 <- ndata$Stime[p1]
        for (k in 1: length(y)) {
          p2 <- which(t_2 <= y[k])
          res.ci[k,j] <- 1 - sum(G0[p2])
          if (res.ci[k,j] < 0) res.ci[k,j] <- 0
        }
        return(res.ci)
      }

      if (isTRUE(cluster)) {
        if (is.null(ncores)) {
          num_cores <- detectCores() - 1
        }else{
          num_cores <- ncores
        }
        registerDoParallel(cores = num_cores)
        on.exit(stopImplicitCluster())

        suppressMessages(
          res.ci <- foreach(i = 1:n.boot, .combine = cbind) %dorng%
            simplebootsurvLDM(object, y, x, lower.tail)
        )

      }else{
        suppressMessages(
          res.ci <- foreach(i = 1:n.boot, .combine = cbind) %do%
            simplebootsurvLDM(object, y, x, lower.tail)
        )
      }


      for (k in 1: length(y)) {
        res.li[k] <- quantile(res.ci[k,], (1 - conf.level) / 2)
        res.ls[k] <- quantile(res.ci[k,], 1 - (1 - conf.level) / 2)
        if (res.li[k] < 0) res.li[k] <- 0
        if (res.ls[k] < 0) res.ls[k] <- 0
      }

    }

    if(conf==FALSE) {
      result <- list(est=resu, estimate=res[!ii], y=y, x=x, conf=conf)
    }
    if(conf==TRUE) {
      result <- list(est=resu, estimate=res[!ii], LCI=res.li, UCI=res.ls, conf.level=conf.level, y=y, x=x, conf=conf)
    }
    class(result) <- c("LDM", "survCS")
    return(invisible(result))
  }
