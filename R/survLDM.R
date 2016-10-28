survLDM <-
  function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95,
           lower.tail = FALSE, cluster = FALSE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    if (missing(x))
      x <- 0
    if (missing(y))
      y <- max(object[[1]]$Stime)
    if (length(x) != length(lower.tail))
      stop("Arguments 'x' and 'lower.tail' must have the same length")
    lenc <- dim(object[[1]])[2]
    ntimes <- lenc%/%2
    if (any(x < 0))
      stop("'x' values should be nonnegative")
    if (any(y < 0))
      stop("'y' values should be nonnegative")
    if (any(y < max(x)))
      stop("'y' values should be equal or greater than all values in 'x'")
    if (length(x) != ntimes-1) {
      cat("The number of consecutive event times in", substitute(object), "is", ntimes, ". The length of 'x' should be", ntimes-1,"\n")
      stop("The length of 'x' is not supported for the selected 'object'")
    }

    y <- y[y >= max(x)]
    y <- sort(unique(y))

    text1 <- paste0("T", c(1:length(x)))
    text2 <- ifelse(lower.tail == TRUE, "<=", ">")
    text3 <- paste0(text1, text2, x, collapse = ",")

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
    res.li <- rep(0, length(y))
    res.ls <- rep(0, length(y))
    resu <- data.frame(cbind(y, res))
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

      if (length(y) == 1) cat("P(T>",y,"|", text3, ") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep = "", "\n")

      if (length(y)>1) {
        resu <- data.frame(cbind(resu,res.li,res.ls))
        names(resu) <- c("y","estimate","LCI","UCI")
        cat("Estimates of ", sep="")
        cat("P(T>y|", text3, ")", sep = "", "\n")
        print(resu)
      }

    }

    if(conf==FALSE) {
      result <- list(est=resu, estimate=res, y=y, x=x, conf=conf)

      if (length(y) == 1) cat("P(T>",y,"|",text3, ") = ", res, sep="", "\n")

      if (length(y)>1) {
        cat("Estimates of ", sep="")
        cat("P(T>y|", text3, ")", sep = "", "\n")
        print(resu)
      }
    }
    if(conf==TRUE) {
      result <- list(est=resu, estimate=res, LCI=res.li, UCI=res.ls, conf.level=conf.level, y=y, x=x, conf=conf)
    }
    class(result) <- c("LDM", "surv")
    return(invisible(result))
  }
