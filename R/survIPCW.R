survIPCW <-
  function(object, x, y, z.name, z.value, bw = "dpik", window = "gaussian",
           method.weights = "NW", conf = FALSE, n.boot = 200,
           conf.level = 0.95, lower.tail = FALSE, cluster = TRUE,
           ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    if (class(object)[1]!="survCS") stop("The argumment 'object' must be of classe 'survCS'")
    obj <- object[[1]]
    if (missing(x))
      x <- 0
    if (missing(y))
      y <- max(object[[1]]$Stime)
    if (missing(z.name) | ncol(object[[1]]) < 5) stop("Argument 'z.name' is missing, with no default")
    covar <- which(names(object[[1]]) == z.name)
    if (missing(z.value)) z.value <- mean(object[[1]][, covar])
    if (length(x) > 1) stop("Length of 'x' must be 1")
    res <- rep(0, length(z.value))
    ifelse (is.numeric(bw), lbd2 <- bw, lbd2 <- dpik(x = object[[1]][, covar]))
    n <- dim(object[[1]])[1]
    delta4 <- rep(1, n)

    y <- y[ y>=x ]
    y <- sort(unique(y))
    z.value <- sort(unique(z.value))
    if (length(y) > 1 & length(z.value) > 1) stop("Argumments 'y' and 'z' cannot have both length greater than 1")

    n.auxi <- max(length(y), length(z.value))
    res <- rep(0, n.auxi)
    res.li <- rep(0, n.auxi)
    res.ls <- rep(0, n.auxi)

    if (length(y) == length(z.value)){
      auxi <- 1   }
    if (length(y) > length(z.value)) {
      auxi <- 2   }
    if (length(y) < length(z.value)) {
      auxi <- 3   }


    if (auxi == 1 | auxi == 3) {
      for (j in 1:n.auxi) {
        be <- rep(0, n)
        be2 <- rep(0, n)
        for (k in 1:n){
          be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event, object[[1]][, covar], delta4, z.value[j], object[[1]]$Stime[k], kernel = window, bw = lbd2)
          be2[k] <- Beran(object[[1]]$time1, 1 - object[[1]]$event1, object[[1]][, covar], delta4, z.value[j], object[[1]]$time1[k], kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(object[[1]][, covar], z.value[j], kernel = window, bw = lbd2), w1 <- LLW(object[[1]][, covar], bw = lbd2, t1 = z.value[j]))
        if (lower.tail == FALSE) {
          p1 <- which(object[[1]]$Stime <= y & object[[1]]$time1 > x & object[[1]]$event == 1)
          p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 == 1)
          ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res[j] <- NA, res[j] <- 1 - sum(w1[p1] / be[p1]) / (1 - sum(w1[p2] / be2[p2])))
        }
        if (lower.tail == TRUE) {
          p1 <- which(object[[1]]$Stime <= y & object[[1]]$time1 <= x & object[[1]]$event == 1)
          p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 == 1)
          ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res[j] <- NA, res[j] <- 1 - sum(w1[p1] / be[p1]) / (sum(w1[p2] / be2[p2])))
        }
        if (sum(w1[p2]) == 0) res[j] <- NA
        if (is.na(res[j])) count <- count + 1
        else {
          if (res[j] > 1) res[j] <- 1
          if (res[j] < 0) res[j] <- 0}
      }

      resu <- data.frame(cbind(z.value, res))
      names(resu) <- c("z", "estimate")
    }

    if (auxi == 2) {
      for (j in 1:n.auxi) {
        be <- rep(0, n)
        be2 <- rep(0, n)
        for (k in 1:n){
          be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event, object[[1]][, covar], delta4, z.value, object[[1]]$Stime[k], kernel = window, bw = lbd2)
          be2[k] <- Beran(object[[1]]$time1, 1 - object[[1]]$event1, object[[1]][, covar], delta4, z.value, object[[1]]$time1[k], kernel = window, bw = lbd2)}
        ifelse(method.weights == "NW", w1 <- NWW(object[[1]][, covar], z.value, kernel = window, bw = lbd2), w1 <- LLW(object[[1]][, covar], bw = lbd2, t1 = z.value))
        if (lower.tail == FALSE) {
          p1 <- which(object[[1]]$Stime <= y[j] & object[[1]]$time1 > x & object[[1]]$event == 1)
          p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 == 1)
          ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res[j] <- NA, res[j] <- 1 - sum(w1[p1] / be[p1]) / (1 - sum(w1[p2] / be2[p2])))
        }
        if (lower.tail == TRUE) {
          p1 <- which(object[[1]]$Stime <= y[j] & object[[1]]$time1 <= x & object[[1]]$event == 1)
          p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 == 1)
          ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res[j] <- NA, res[j] <- 1 - sum(w1[p1] / be[p1]) / (sum(w1[p2] / be2[p2])))
        }
        if (sum(w1[p2]) == 0) res[j] <- NA
        if (is.na(res[j])) count <- count + 1
        else {
          if (res[j] > 1) res[j] <- 1
          if (res[j] < 0) res[j] <- 0}
      }
      resu <- data.frame(cbind(y, res))
      names(resu) <- c("y", "estimate")
    }

    if (conf == TRUE) {

      simpleboot <- function(object, n, n.auxi, auxi){
        k <- 1
        res.ci <- matrix(0, nrow = n.auxi, ncol = k)
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        obj <- ndata
        ifelse (is.numeric(bw), lbd2 <- bw, lbd2 <- dpik(x = ndata[, covar]))

        if (auxi == 1 | auxi == 3) {
          for (j in 1:n.auxi) {
            be <- rep(0, n)
            be2 <- rep(0, n)
            for (i in 1:n) {
              be[i] <- Beran(ndata$Stime, 1 - ndata$event, ndata[, covar], delta4, z.value[j], ndata$Stime[i], kernel = window, bw = lbd2)
              be2[i] <- Beran(ndata$time1, 1 - ndata$event1, ndata[, covar], delta4, z.value[j], ndata$time1[i], kernel = window, bw = lbd2)
            }
            ifelse(method.weights == "NW",
                   w1 <- NWW(ndata[, covar], z.value[j], kernel = window, bw = lbd2),
                   w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value[j]))
            if (lower.tail == FALSE) {
              p1 <- which(ndata$Stime <= y & ndata$time1 > x & ndata$event == 1)
              p2 <- which(ndata$time1 <= x & ndata$event1 == 1)
              ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 - sum(w1[p1] / be[p1]) / (1 - sum(w1[p2] / be2[p2])))
            }
            if (lower.tail == TRUE) {
              p1 <- which(ndata$Stime <= y & ndata$time1 <= x & ndata$event == 1)
              p2 <- which(ndata$time1 <= x & ndata$event1 == 1)
              ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 - sum(w1[p1] / be[p1]) / (sum(w1[p2] / be2[p2])))
            }
            if (sum(w1[p2]) == 0) res.ci[j, k] <- NA
            if (is.na(res.ci[j, k])) count <- count + 1
            else {
              if (res.ci[j, k] > 1) res.ci[j, k] <- 1
              if (res.ci[j, k] < 0) res.ci[j, k] <- 0}
          }
        }

        if (auxi == 2) {
          for (j in 1:n.auxi) {
            be <- rep(0, n)
            be2 <- rep(0, n)
            for (i in 1:n){
              be[i] <- Beran(ndata$Stime, 1 - ndata$event, ndata[, covar], delta4, z.value, ndata$Stime[i], kernel = window, bw = lbd2)
              be2[i] <- Beran(ndata$time1, 1 - ndata$event1, ndata[, covar], delta4, z.value, ndata$time1[i], kernel = window, bw = lbd2)
            }
            ifelse(method.weights == "NW", w1 <- NWW(ndata[, covar], z.value, kernel = window, bw = lbd2), w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value))
            if (lower.tail == FALSE) {
              p1 <- which(ndata$Stime <= y[j] & ndata$time1 > x & ndata$event == 1)
              p2 <- which(ndata$time1 <= x & ndata$event1 == 1)
              ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 - sum(w1[p1] / be[p1]) / (1 - sum(w1[p2] / be2[p2])))
            }
            if (lower.tail == TRUE) {
              p1 <- which(ndata$Stime <= y[j] & ndata$time1 <= x & ndata$event == 1)
              p2 <- which(ndata$time1 <= x & ndata$event1 == 1)
              ifelse (any(be[p1] == 0) | any(be2[p2] == 0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 - sum(w1[p1] / be[p1]) / (sum(w1[p2] / be2[p2])))
            }
            if (sum(w1[p2]) == 0) res.ci[j, k] <- NA
            if (is.na(res.ci[j, k])) count <- count + 1
            else {
              if (res.ci[j, k] > 1) res.ci[j, k] <- 1
              if (res.ci[j, k] < 0) res.ci[j, k] <- 0}
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
            simpleboot(object, n, n.auxi, auxi)
        )

      }else{
        res.ci <- foreach(i = 1:n.boot, .combine = cbind,
                          .export = "simpleboot") %do%
          simpleboot(object, n, n.auxi, auxi)

      }



      for (k in 1: n.auxi) {
        res.li[k] <- quantile(res.ci[k,], (1 - conf.level) / 2, na.rm = TRUE)
        res.ls[k] <- quantile(res.ci[k,], 1 - ( 1 - conf.level) / 2, na.rm = TRUE)
      }

      if (auxi == 1 | auxi == 3) { resu <- data.frame(cbind(z.value, res, res.li, res.ls))
      names(resu) <- c("z", "estimate", "LCI", "UCI")}
      if (auxi == 2) { resu <- data.frame(cbind(y, res, res.li, res.ls))
      names(resu) <- c("y", "estimate", "LCI", "UCI")}
    }

    if(conf == FALSE & auxi == 1) {
      result <- list(est = resu, estimate = res, z.name = z.name, z.value = z.value, x = x, y = y, bw = bw, window = window, method.weights = method.weights, conf = conf, lbd = lbd2)
      cat("S(T>",y,"|T1>",x,",",z.name,"=",z.value,") = ", res, sep="", "\n")
    }

    if(conf == FALSE & auxi == 3) {
      result <- list(est = resu, estimate = res, z.name = z.name, z.value = z.value, x = x, y = y, bw = bw, window = window, method.weights = method.weights, conf = conf, lbd = lbd2)
      cat("Estimates of S(T>",y,"|T1>",x,",",z.name,")", sep="", "\n")
      print(resu)
    }

    if(conf == FALSE & auxi == 2) {
      result <- list(est = resu, estimate = res, z.name = z.name, z.value = z.value, x = x, y = y, bw = bw, window = window, method.weights = method.weights, conf = conf, lbd = lbd2)
      cat("Estimates of S(T> y|T1>",x,",",z.name,"=",z.value,")", sep="", "\n")
      print(resu)
    }

    if(conf == TRUE & auxi == 1) {
      result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, z.name = z.name, z.value = z.value, bw = bw, window = window, method.weights = method.weights, x = x, y = y, conf = conf, lbd = lbd2)
      cat("S(T>",y,"|T1>",x,",",z.name,"=",z.value,") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
    }

    if(conf == TRUE & auxi == 2) {
      result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, z.name = z.name, z.value = z.value, bw = bw, window = window, method.weights = method.weights, x = x, y = y, conf = conf, lbd = lbd2)
      cat("Estimates of S(T> y|T1>",x,",",z.name,"=",z.value,")", sep="", "\n")
      print(resu)
    }

    if(conf == TRUE & auxi == 3) {
      result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, z.name = z.name, z.value = z.value, bw = bw, window = window, method.weights = method.weights, conf = conf, lbd = lbd2)
      cat("Estimates of S(T>",y,"|T1>",x,",",z.name,")", sep="", "\n")
      print(resu)
    }

    class(result) <- "IPCW"
    return(invisible(result))
  }
