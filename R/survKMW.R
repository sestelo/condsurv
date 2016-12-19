survKMW <-
  function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95,
           lower.tail = FALSE, cluster = FALSE, ncores = NULL, na.rm=T)
  {


    lenc <- dim(object[[1]])[2]
    ntimes <- lenc%/%2

    res <- rep(0, length(y))


    if (ntimes == 2) {
      if (lower.tail == FALSE)
      {
        if(x == 0) {
          for (k in 1: length(y)) { res[k] <- KM(object[[1]]$Stime, object[[1]]$event, y[k])}
        }
        else{
          G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
          G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
          p2 <- which(object[[1]]$time1 <= x)
          if (length(p2) == 0){
            p2 <- which(object[[1]]$time1 > x)
            den <- 1
          }
          else { den <- 1 - sum(G2[p2])}

          if (den == 0) stop("Insufficient data.")
          for (k in 1: length(y)) {
            p1 <- which(object[[1]]$time1 > x & object[[1]]$Stime <= y[k])
            if (length(p1) == 0) res[k] <- 1
            else res[k] <- 1 - sum(G1[p1])/den
            if (res[k] < 0) res[k] <- 0
          }
        }
      }

      if (lower.tail == TRUE)
      {
        G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
        G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
        p2 <- which(object[[1]]$time1 <= x)
        if (length(p2) == 0) stop("Insufficient data.")
        for (k in 1: length(y)) {
          p1 <- which(object[[1]]$time1 <= x & object[[1]]$Stime <= y[k])
          if (length(p1) == 0) { num <- 0
          }else {num <- sum(G1[p1])}
          den <- sum(G2[p2])
          if (den == 0){
            res[k] <- NA
          }else{
            res[k] <- 1 - num/ den
            if (res[k] < 0) res[k] <- 0
          }
        }
      }
    }

    if (ntimes > 2) {

      G2 <- KMW(object[[1]][,2*ntimes-3], object[[1]][,2*ntimes-2])
      G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
      X <- data.frame(object[[1]][,2*(1:ntimes)-1])
      p2 <- whichCS(X, x=x, lower.tail=lower.tail)
      if (length(p2) == 0) stop("Insufficient data.")
      for (k in 1: length(y)) {
        X <- data.frame(object[[1]][,2*(1:ntimes)-1])
        xy <- c(x,y[k])
        lower.tail.y <- c(lower.tail,TRUE)
        p1 <- whichCS(X, x=xy, lower.tail=lower.tail.y)
        if (length(p1) == 0) { num <- 0
        } else {num <- sum(G1[p1])}
        den <- sum(G2[p2])
        if (den == 0){
          res[k] <- NA
        }else{
          res[k] <- 1 - num/ den
          if (res[k] < 0) res[k] <- 0
        }
      }

    }

    ii <- duplicated(res)
    y <- y[!ii]

    res.li <- rep(0, length(y))
    res.ls <- rep(0, length(y))
    resu <- data.frame(cbind(y, res[!ii]))
    names(resu) <- c("y", "estimate")

    if (conf == TRUE) {

      simplebootsurvKMW <- function(object, lower.tail = lower.tail, y, x){
        j <- 1
        res.ci <- matrix(0, nrow = length(y), ncol = j)
        n <- dim(object[[1]])[1]
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]

        if (ntimes == 2) {
          if (lower.tail == FALSE)
          {
            if(x == 0) {
              for (k in 1: length(y)) { res.ci[k,j] <- KM(ndata$Stime, ndata$event, y[k])}
            }
            else{
              G2 <- KMW(ndata$time1, ndata$event1)
              G1 <- KMW(ndata$Stime, ndata$event)
              p2 <- which(ndata$time1 <= x)

              if (length(p2) == 0){
                p2 <- which(ndata$time1 > x)
                den <- 1
              }
              else { den <- 1 - sum(G2[p2])}

              if (den == 0) stop("Insufficient data.")


              #if (length(p2) == 0) stop("Insufficient data.")
              for (k in 1: length(y)) {
                p1 <- which(ndata$time1 > x & ndata$Stime <= y[k])
                if (length(p1) == 0) res.ci[k, j] <- 1
                else res.ci[k, j] <- 1 - sum(G1[p1]) / den
                if (res.ci[k, j] < 0) res.ci[k, j] <- 0
              }
            }
          }

          if (lower.tail == TRUE)
          {
            G2 <- KMW(ndata$time1, ndata$event1)
            G1 <- KMW(ndata$Stime, ndata$event)
            p2 <- which(ndata$time1 <= x)
            if (length(p2) == 0) stop("Insufficient data.")
            for (k in 1: length(y)) {
              p1 <- which(ndata$time1 <= x & ndata$Stime <= y[k])
              if (length(p1) == 0){ num <- 0
              }else{ num <- sum(G1[p1])}
              den <- sum(G2[p2])
              if (den == 0){
                res.ci[k, j] <- NA
              }else{
                res.ci[k, j] <- 1 - num / den
                if (res.ci[k, j] < 0) res.ci[k,j] <- 0
              }
            }
          }
        }

        if (ntimes > 2) {
          G2 <- KMW(ndata[,2*ntimes-3], ndata[,2*ntimes-2])
          G1 <- KMW(ndata$Stime, ndata$event)
          #p2 <- which(ndata$time1 <= x)
          X <- data.frame(ndata[,2*(1:ntimes)-1])
          p2 <- whichCS(X, x=x, lower.tail=lower.tail)
          if (length(p2) == 0) stop("Insufficient data.")

          for (k in 1: length(y)) {
            #p1 <- which(ndata$time1 <= x & ndata$Stime <= y[k])
            X <- data.frame(object[[1]][,2*(1:ntimes)-1])
            xy <- c(x,y[k])
            lower.tail.y <- c(lower.tail,TRUE)
            p1 <- whichCS(X, x=xy, lower.tail=lower.tail.y)
            if (length(p1) == 0) { num <- 0
            } else {num <- sum(G1[p1])}
            den <- sum(G2[p2])
            if (den == 0){
              res.ci[k, j] <- NA
            }else{
              res.ci[k, j] <- 1 - num / den
              if (res.ci[k, j] < 0) res.ci[k,j] <- 0
            }
          }
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
            simplebootsurvKMW(object, lower.tail, y, x)
        )

      }else{
        suppressMessages(
          res.ci <- foreach(i = 1:n.boot, .combine = cbind) %do%
            simplebootsurvKMW(object, lower.tail, y, x)
        )
      }

      for (k in 1: length(y)) {
        res.li[k] <- quantile(res.ci[k,], (1 - conf.level) / 2, na.rm=na.rm)
        res.ls[k] <- quantile(res.ci[k,], 1 - (1 - conf.level) / 2, na.rm=na.rm)
      }

      if (length(y)>1) {
        resu <- data.frame(cbind(resu,res.li,res.ls))
        names(resu) <- c("y","estimate","LCI","UCI")
      }

    }

    if(conf == FALSE) { result <- list(est=resu, estimate=res[!ii], y=y, x=x, conf=conf) }

    if(conf == TRUE) { result <- list(est = resu, estimate = res[!ii], LCI = res.li, UCI = res.ls, conf.level = conf.level, y = y, x = x, conf = conf) }

    class(result) <- c("KMW", "survCS")
    return(invisible(result))
  }

