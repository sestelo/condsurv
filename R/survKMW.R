survKMW <-
  function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95,
           lower.tail = FALSE, cluster = TRUE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    if (missing(x))
      x <- 0
    if (missing(y))
      y <- max(object[[1]]$Stime)
    if (length(x) > 1) stop("Length of 'x' must be 1")

    y <- y[y >= x]
    y <- sort(unique(y))
    res <- rep(0, length(y))
    res.li <- rep(0, length(y))
    res.ls <- rep(0, length(y))


    if (lower.tail == FALSE)
    {
      G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
      G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
      p2 <- which(object[[1]]$time1 <= x)
      for (k in 1: length(y)) {
        p1 <- which(object[[1]]$time1 > x & object[[1]]$Stime <= y[k])
        res[k] <- 1 - sum(G1[p1])/(1 - sum(G2[p2]))
      }
    }

    if (lower.tail == TRUE)
    {
      G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
      G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
      p2 <- which(object[[1]]$time1 <= x)
      for (k in 1: length(y)) {
        p1 <- which(object[[1]]$time1 <= x & object[[1]]$Stime <= y[k])
        res[k] <- 1 - sum(G1[p1])/ sum(G2[p2])
      }
    }

    resu <- data.frame(cbind(y, res))
    names(resu) <- c("y", "estimate")

    if (conf == TRUE) {

      simpleboot <- function(object, lower.tail = lower.tail, y, x){
        j <- 1
        res.ci <- matrix(0, nrow = length(y), ncol = j)
        n <- dim(object[[1]])[1]
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]


        if (lower.tail == FALSE)
        {
          G2 <- KMW(ndata$time1, ndata$event1)
          G1 <- KMW(ndata$Stime, ndata$event)
          p2 <- which(ndata$time1 <= x)
          for (k in 1: length(y)) {
            p1 <- which(ndata$time1 > x & ndata$Stime <= y[k])
            res.ci[k, j] <- 1 - sum(G1[p1]) / (1 - sum(G2[p2]))
          }
        }

        if (lower.tail == TRUE)
        {
          G2 <- KMW(ndata$time1, ndata$event1)
          G1 <- KMW(ndata$Stime, ndata$event)
          p2 <- which(ndata$time1 <= x)
          for (k in 1: length(y)) {
            p1 <- which(ndata$time1 <= x & ndata$Stime <= y[k])
            res.ci[k, j] <- 1 - sum(G1[p1]) / sum(G2[p2])
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
          res.ci <- foreach(i = 1:n.boot, .combine = cbind,
                            .export = "simpleboot") %dorng%
            simpleboot(object, lower.tail, y, x)
        )

      }else{
        res.ci <- foreach(i = 1:n.boot, .combine = cbind,
                          .export = "simpleboot") %do%
          simpleboot(object, lower.tail, y, x)

      }

      for (k in 1: length(y)) {
        res.li[k] <- quantile(res.ci[k,], (1 - conf.level) / 2)
        res.ls[k] <- quantile(res.ci[k,], 1 - (1 - conf.level) / 2)
      }
      if (length(y) == 1 & lower.tail == FALSE) cat("S(T>",y,"|T1>",x,") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
      if (length(y) == 1 & lower.tail == TRUE) cat("S(T>",y,"|T1<=",x,") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
      if (length(y) > 1) {
        resu <- data.frame(cbind(resu, res.li, res.ls))
        names(resu) <- c("y", "estimate", "LCI", "UCI")
        if (lower.tail == FALSE) cat("Estimates of S(T>y|T1>",x,")", sep = "", "\n")
        if (lower.tail == TRUE) cat("Estimates of S(T>y|T1<=",x,")", sep = "", "\n")
        print(resu)
      }
    }

    if(conf==FALSE) {
      result <- list(est = resu, estimate = res, y = y, x = x, conf = conf)
      if (length(y) == 1 & lower.tail == FALSE) cat("S(T>",y,"|T1>",x,") = ", res, sep="", "\n")
      if (length(y) == 1 & lower.tail == TRUE) cat("S(T>",y,"|T1<=",x,") = ", res, sep="", "\n")
      if (length(y) > 1) {
        if (lower.tail == FALSE) cat("Estimates of S(T>y|T1>",x,")",sep="","\n")
        if (lower.tail == TRUE) cat("Estimates of S(T>y|T1<=",x,")",sep="","\n")
        print(resu)
      }
    }

    if(conf == TRUE) { result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, y = y, x = x, conf = conf) }

    class(result) <- c("KMW", "surv")
    return(invisible(result))
  }
