survKMW <-
  function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95,
           lower.tail = FALSE, cluster = FALSE, ncores = NULL, na.rm=T)
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
    
    y <- y[y >= max(x)]
    y <- sort(unique(y))
    
    text1 <- paste("T",c(1:length(x)),sep="")
    text2 <- ifelse(lower.tail==TRUE, "<=",">")
    text3 <- paste(text1,text2,x, sep="")
    
    res <- rep(0, length(y))
    res.li <- rep(0, length(y))
    res.ls <- rep(0, length(y))
    
    if (ntimes == 2) {
      if (lower.tail == FALSE)
      {
        G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
        G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
        p2 <- which(object[[1]]$time1 <= x)
        for (k in 1: length(y)) {
          p1 <- which(object[[1]]$time1 > x & object[[1]]$Stime <= y[k])
          res[k] <- 1 - sum(G1[p1])/(1 - sum(G2[p2]))
          if (res[k] < 0) res[k] <- 0
        }
      }
      
      if (lower.tail == TRUE)
      {
        G2 <- KMW(object[[1]]$time1, object[[1]]$event1)
        G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
        p2 <- which(object[[1]]$time1 <= x)
        for (k in 1: length(y)) {
          p1 <- which(object[[1]]$time1 <= x & object[[1]]$Stime <= y[k])
          den <- sum(G2[p2])
          res[k] <- 1 - sum(G1[p1])/ den
          if (res[k] < 0) res[k] <- 0
          if (den == 0) res[k] <- NA
        }
      }
    }
    
    if (ntimes > 2) {
      
      G2 <- KMW(object[[1]][,2*ntimes-3], object[[1]][,2*ntimes-2])
      G1 <- KMW(object[[1]]$Stime, object[[1]]$event)
      #p2 <- which(object[[1]]$time1 <= x)
      X <- data.frame(object[[1]][,2*(1:ntimes)-1])
      p2 <- whichCS(X, x=x, lower.tail=lower.tail)
      
      for (k in 1: length(y)) {
        #p1 <- which(object[[1]]$time1 <= x & object[[1]]$Stime <= y[k])
        X <- data.frame(object[[1]][,2*(1:ntimes)-1])
        xy <- c(x,y[k])
        lower.tail.y <- c(lower.tail,TRUE)
        p1 <- whichCS(X, x=xy, lower.tail=lower.tail.y)
        den <- sum(G2[p2])
        res[k] <- 1 - sum(G1[p1])/ den
        if (res[k] < 0) res[k] <- 0
        if (den == 0) res[k] <- NA
      }
      
    }
    
    resu <- data.frame(cbind(y, res))
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
            G2 <- KMW(ndata$time1, ndata$event1)
            G1 <- KMW(ndata$Stime, ndata$event)
            p2 <- which(ndata$time1 <= x)
            for (k in 1: length(y)) {
              p1 <- which(ndata$time1 > x & ndata$Stime <= y[k])
              res.ci[k, j] <- 1 - sum(G1[p1]) / (1 - sum(G2[p2]))
              if (res.ci[k, j] < 0) res.ci[k, j] <- 0
            }
          }
          
          if (lower.tail == TRUE)
          {
            G2 <- KMW(ndata$time1, ndata$event1)
            G1 <- KMW(ndata$Stime, ndata$event)
            p2 <- which(ndata$time1 <= x)
            for (k in 1: length(y)) {
              p1 <- which(ndata$time1 <= x & ndata$Stime <= y[k])
              den <- sum(G2[p2]) 
              res.ci[k, j] <- 1 - sum(G1[p1]) / den
              if (res.ci[k, j] < 0) res.ci[k,j] <- 0
              if (den == 0) res.ci[k, j] <- NA
            }
          }
        }
        
        if (ntimes > 2) {
          G2 <- KMW(ndata[,2*ntimes-3], ndata[,2*ntimes-2])
          G1 <- KMW(ndata$Stime, ndata$event)
          #p2 <- which(ndata$time1 <= x)
          X <- data.frame(ndata[,2*(1:ntimes)-1])
          p2 <- whichCS(X, x=x, lower.tail=lower.tail)
          
          for (k in 1: length(y)) {
            #p1 <- which(ndata$time1 <= x & ndata$Stime <= y[k])
            X <- data.frame(object[[1]][,2*(1:ntimes)-1])
            xy <- c(x,y[k])
            lower.tail.y <- c(lower.tail,TRUE)
            p1 <- whichCS(X, x=xy, lower.tail=lower.tail.y)
            den <- sum(G2[p2])
            res.ci[k, j] <- 1 - sum(G1[p1]) / den
            if (res.ci[k, j] < 0) res.ci[k,j] <- 0
            if (den == 0) res.ci[k, j] <- NA
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
      
      if (length(y) == 1) cat(cat("P(T>",y,"|",sep=""), cat(text3, sep=","),") =",res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
      if (length(y)>1) {
        resu <- data.frame(cbind(resu, res.li, res.ls))
        names(resu) <- c("y", "estimate", "LCI", "UCI")
        cat("Estimates of ", sep="")
        cat(cat("P(T>y|",sep=""), cat(text3, sep=","),")",sep="","\n")
        print(resu)
      }
    }
    
    if(conf==FALSE) {
      result <- list(est = resu, estimate = res, y = y, x = x, conf = conf)
      if (length(y) == 1) cat(cat("P(T>",y,"|",sep=""), cat(text3, sep=","),") =",res, sep="", "\n")
      
      if (length(y)>1) {
        cat("Estimates of ", sep="")
        cat(cat("P(T>y|",sep=""), cat(text3, sep=","),")",sep="","\n")
        print(resu)
      }
    }
    
    if(conf == TRUE) { result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, y = y, x = x, conf = conf) }
    
    class(result) <- c("KMW", "surv")
    return(invisible(result))
  }