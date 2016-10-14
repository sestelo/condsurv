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
    
    y <- y[y >= max(x)]
    y <- sort(unique(y))
    
    text1 <- paste("T",c(1:length(x)),sep="")
    text2 <- ifelse(lower.tail==TRUE, "<=",">")
    text3 <- paste(text1,text2,x, sep="")
    
    X <- data.frame(object[[1]][,2*(1:ntimes)-1])
    p1 <- whichCS(X, x=x, lower.tail=lower.tail)
    G0 <- KMW(object[[1]]$Stime[p1], object[[1]]$event[p1])
    t_2 <- object[[1]]$Stime[p1]
    res <- rep(0, length(y))
    for (k in 1: length(y)) {
      p2 <- which(t_2 <= y[k])
      res[k] <- 1 - sum(G0[p2])
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
        G0 <- KMW(ndata$Stime[p1], ndata$event[p1])
        t_2 <- ndata$Stime[p1]
        for (k in 1: length(y)) {
          p2 <- which(t_2 <= y[k])
          res.ci[k,j] <- 1 - sum(G0[p2])
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
      }
      
      if (length(y) == 1) cat(cat("P(T>",y,"|",sep=""), cat(text3, sep=","),") =",res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
      
      if (length(y)>1) {
        resu <- data.frame(cbind(resu,res.li,res.ls))
        names(resu) <- c("y","estimate","LCI","UCI")
        cat("Estimates of ", sep="")
        cat(cat("P(T>y|",sep=""), cat(text3, sep=","),")",sep="","\n")
        print(resu)
      }
      
    }
    
    if(conf==FALSE) {
      result <- list(est=resu, estimate=res, y=y, x=x, conf=conf)
      
      if (length(y) == 1) cat(cat("P(T>",y,"|",sep=""), cat(text3, sep=","),") =",res, sep="", "\n")
      
      if (length(y)>1) {
        cat("Estimates of ", sep="")
        cat(cat("P(T>y|",sep=""), cat(text3, sep=","),")",sep="","\n")
        print(resu)
      }
    }
    if(conf==TRUE) {
      result <- list(est=resu, estimate=res, LCI=res.li, UCI=res.ls, conf.level=conf.level, y=y, x=x, conf=conf)
    }
    class(result) <- c("LDM", "surv")
    return(invisible(result))
  }
