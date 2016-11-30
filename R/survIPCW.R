survIPCW <-
  function (object, x, y, z.name, z.value, bw = "dpik", window = "gaussian",
            method.weights = "NW", conf = FALSE, n.boot = 200, conf.level = 0.95,
            lower.tail = FALSE, cluster = FALSE, ncores = NULL)
  {


    obj <- object[[1]]

    covar <- which(names(object[[1]]) == z.name)
    if (missing(z.value))
      z.value <- mean(object[[1]][, covar])
    res <- rep(0, length(z.value))
    if (bw == "dpik") {
      lbd2 <- dpik(x = object[[1]][, covar])
    }
    else if (bw == "np") {
      options(np.messages = FALSE)
      lbd2 <- npudensbw(dat = object[[1]][, covar])$bw
    }
    else {
      lbd2 <- bw
    }
    if (!is.numeric(bw) & !(bw %in% c("dpik", "np"))) {
      stop("Argument 'bw' have to be 'dpik', 'np' or a numeric.")
    }
    n <- dim(object[[1]])[1]
    delta4 <- rep(1, n)
    lenc <- dim(object[[1]])[2]
    ntimes <- lenc%/%2


    if(length(z.value) > 1){
      stop("Argumment 'z' cannot have length greater than 1")
    }



    n.auxi <- length(y)

    res <- rep(0, n.auxi)



    if (length(y) == length(z.value)) {
      auxi <- 1
    }
    if (length(y) > length(z.value)) {
      auxi <- 2
    }




    if (ntimes == 2) {
      if (auxi == 1) {
        for (j in 1:n.auxi) {
          be <- rep(0, n)
          be2 <- rep(0, n)
          for (k in 1:n) {
            be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event,
                           object[[1]][, covar], delta4, z.value[j],
                           object[[1]]$Stime[k], kernel = window, bw = lbd2)
            be2[k] <- Beran(object[[1]]$time1, 1 - object[[1]]$event1,
                            object[[1]][, covar], delta4, z.value[j],
                            object[[1]]$time1[k], kernel = window, bw = lbd2)
          }
          ifelse(method.weights == "NW", w1 <- NWW(object[[1]][,
                                                               covar], z.value[j], kernel = window, bw = lbd2),
                 w1 <- LLW(object[[1]][, covar], bw = lbd2,
                           t1 = z.value[j]))

          if (lower.tail == FALSE) {
            p1 <- which(object[[1]]$Stime <= y & object[[1]]$time1 >
                          x & object[[1]]$event == 1)
            p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 ==
                          1)
            ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                   res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(1 -
                                                                     sum(w1[p2]/be2[p2])))
          }

          if (lower.tail == TRUE) {
            p1 <- which(object[[1]]$Stime <= y & object[[1]]$time1 <=
                          x & object[[1]]$event == 1)
            p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 ==
                          1)
            ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                   res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
          }

          if (sum(w1[p2]) == 0)
            res[j] <- NA
          if (is.na(res[j]))
            count <- count + 1
          else {
            if (res[j] > 1)
              res[j] <- 1
            if (res[j] < 0)
              res[j] <- 0
          }
        }
        resu <- data.frame(cbind(y, res))
        names(resu) <- c("y", "estimate")
      }

      if (auxi == 2) {

        be <- rep(0, n)
        be2 <- rep(0, n)
        for (k in 1:n) {
          be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event,
                         object[[1]][, covar], delta4, z.value, object[[1]]$Stime[k],
                         kernel = window, bw = lbd2)
          be2[k] <- Beran(object[[1]]$time1, 1 - object[[1]]$event1,
                          object[[1]][, covar], delta4, z.value, object[[1]]$time1[k],
                          kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(object[[1]][,
                                                             covar], z.value, kernel = window, bw = lbd2),
               w1 <- LLW(object[[1]][, covar], bw = lbd2, t1 = z.value))

        for (j in 1:n.auxi) {

          if (lower.tail == FALSE) {
            p1 <- which(object[[1]]$Stime <= y[j] & object[[1]]$time1 >
                          x & object[[1]]$event == 1)
            p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 ==
                          1)
            ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                   res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(1 -
                                                                     sum(w1[p2]/be2[p2])))
          }
          if (lower.tail == TRUE) {
            p1 <- which(object[[1]]$Stime <= y[j] & object[[1]]$time1 <=
                          x & object[[1]]$event == 1)
            p2 <- which(object[[1]]$time1 <= x & object[[1]]$event1 ==
                          1)
            ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                   res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
          }
          if (sum(w1[p2]) == 0)
            res[j] <- NA
          if (is.na(res[j]))
            count <- count + 1
          else {
            if (res[j] > 1)
              res[j] <- 1
            if (res[j] < 0)
              res[j] <- 0
          }
        }
        resu <- data.frame(cbind(y, res))
        names(resu) <- c("y", "estimate")
      }
    }
    if (ntimes > 2) {
      if (auxi == 1) {
        for (j in 1:n.auxi) {
          be <- rep(0, n)
          be2 <- rep(0, n)
          for (k in 1:n) {
            be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event,
                           object[[1]][, covar], delta4, z.value[j],
                           object[[1]]$Stime[k], kernel = window, bw = lbd2)
            be2[k] <- Beran(object[[1]][, 2 * ntimes -
                                          3], 1 - object[[1]][, 2 * ntimes - 2], object[[1]][,
                                                                                             covar], delta4, z.value[j], object[[1]][k,
                                                                                                                                     2 * ntimes - 3], kernel = window, bw = lbd2)
          }
          ifelse(method.weights == "NW", w1 <- NWW(object[[1]][,
                                                               covar], z.value[j], kernel = window, bw = lbd2),
                 w1 <- LLW(object[[1]][, covar], bw = lbd2,
                           t1 = z.value[j]))
          X <- data.frame(object[[1]][, 2 * (1:ntimes) -
                                        1])
          pos1 <- whichCS(X, x = x, lower.tail = lower.tail)
          pos2 <- which(object[[1]]$Stime <= y & object[[1]]$event ==
                          1)
          p1 <- intersect(pos1, pos2)
          pos3 <- which(object[[1]][, 2 * ntimes - 2] ==
                          1)
          p2 <- intersect(pos1, pos3)
          ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                 res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
          if (sum(w1[p2]) == 0)
            res[j] <- NA
          if (is.na(res[j]))
            count <- count + 1
          else {
            if (res[j] > 1)
              res[j] <- 1
            if (res[j] < 0)
              res[j] <- 0
          }
        }
        resu <- data.frame(cbind(y, res))
        names(resu) <- c("y", "estimate")
      }
      if (auxi == 2) {

        be <- rep(0, n)
        be2 <- rep(0, n)
        for (k in 1:n) {
          be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event,
                         object[[1]][, covar], delta4, z.value[j],
                         object[[1]]$Stime[k], kernel = window, bw = lbd2)
          be2[k] <- Beran(object[[1]][, 2 * ntimes -
                                        3], 1 - object[[1]][, 2 * ntimes - 2], object[[1]][,
                                                                                           covar], delta4, z.value[j], object[[1]][k,
                                                                                                                                   2 * ntimes - 3], kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(object[[1]][,
                                                             covar], z.value, kernel = window, bw = lbd2),
               w1 <- LLW(object[[1]][, covar], bw = lbd2, t1 = z.value))

        for (j in 1:n.auxi) {
          X <- data.frame(object[[1]][, 2 * (1:ntimes) -
                                        1])
          pos1 <- whichCS(X, x = x, lower.tail = lower.tail)
          pos2 <- which(object[[1]]$Stime <= y & object[[1]]$event ==
                          1)
          p1 <- intersect(pos1, pos2)
          pos3 <- which(object[[1]][, 2 * ntimes - 2] ==
                          1)
          p2 <- intersect(pos1, pos3)
          ifelse(any(be[p1] == 0) | any(be2[p2] == 0),
                 res[j] <- NA, res[j] <- 1 - sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
          if (sum(w1[p2]) == 0)
            res[j] <- NA
          if (is.na(res[j]))
            count <- count + 1
          else {
            if (res[j] > 1)
              res[j] <- 1
            if (res[j] < 0)
              res[j] <- 0
          }
        }
        resu <- data.frame(cbind(y, res))
        names(resu) <- c("y", "estimate")
      }
    }


    ii <- duplicated(resu$estimate)
    y <- y[!ii]

    n.auxi <- length(y)
    res.li <- rep(0, n.auxi)
    res.ls <- rep(0, n.auxi)


    if (conf == TRUE) {
      simplebootsurvIPCW <- function(object, n, n.auxi, auxi) {
        k <- 1
        res.ci <- matrix(0, nrow = n.auxi, ncol = k)
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx, ]
        obj <- ndata
        ifelse(is.numeric(bw), lbd2 <- bw, lbd2 <- dpik(x = ndata[,
                                                                  covar]))
        if (ntimes == 2) {
          if (auxi == 1) {
            for (j in 1:n.auxi) {
              be <- rep(0, n)
              be2 <- rep(0, n)
              for (i in 1:n) {
                be[i] <- Beran(ndata$Stime, 1 - ndata$event,
                               ndata[, covar], delta4, z.value[j], ndata$Stime[i],
                               kernel = window, bw = lbd2)
                be2[i] <- Beran(ndata$time1, 1 - ndata$event1,
                                ndata[, covar], delta4, z.value[j], ndata$time1[i],
                                kernel = window, bw = lbd2)
              }
              ifelse(method.weights == "NW", w1 <- NWW(ndata[,
                                                             covar], z.value[j], kernel = window, bw = lbd2),
                     w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value[j]))
              if (lower.tail == FALSE) {
                p1 <- which(ndata$Stime <= y & ndata$time1 >
                              x & ndata$event == 1)
                p2 <- which(ndata$time1 <= x & ndata$event1 ==
                              1)
                ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                                0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                         sum(w1[p1]/be[p1])/(1 - sum(w1[p2]/be2[p2])))
              }
              if (lower.tail == TRUE) {
                p1 <- which(ndata$Stime <= y & ndata$time1 <=
                              x & ndata$event == 1)
                p2 <- which(ndata$time1 <= x & ndata$event1 ==
                              1)
                ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                                0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                         sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
              }
              if (sum(w1[p2]) == 0)
                res.ci[j, k] <- NA
              if (is.na(res.ci[j, k]))
                count <- count + 1
              else {
                if (res.ci[j, k] > 1)
                  res.ci[j, k] <- 1
                if (res.ci[j, k] < 0)
                  res.ci[j, k] <- 0
              }
            }
          }
          if (auxi == 2) {

            be <- rep(0, n)
            be2 <- rep(0, n)
            for (i in 1:n) {
              be[i] <- Beran(ndata$Stime, 1 - ndata$event,
                             ndata[, covar], delta4, z.value, ndata$Stime[i],
                             kernel = window, bw = lbd2)
              be2[i] <- Beran(ndata$time1, 1 - ndata$event1,
                              ndata[, covar], delta4, z.value, ndata$time1[i],
                              kernel = window, bw = lbd2)
            }
            ifelse(method.weights == "NW", w1 <- NWW(ndata[,
                                                           covar], z.value, kernel = window, bw = lbd2),
                   w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value))

            for (j in 1:n.auxi) {
              if (lower.tail == FALSE) {
                p1 <- which(ndata$Stime <= y[j] & ndata$time1 >
                              x & ndata$event == 1)
                p2 <- which(ndata$time1 <= x & ndata$event1 ==
                              1)
                ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                                0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                         sum(w1[p1]/be[p1])/(1 - sum(w1[p2]/be2[p2])))
              }
              if (lower.tail == TRUE) {
                p1 <- which(ndata$Stime <= y[j] & ndata$time1 <=
                              x & ndata$event == 1)
                p2 <- which(ndata$time1 <= x & ndata$event1 ==
                              1)
                ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                                0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                         sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
              }
              if (sum(w1[p2]) == 0)
                res.ci[j, k] <- NA
              if (is.na(res.ci[j, k]))
                count <- count + 1
              else {
                if (res.ci[j, k] > 1)
                  res.ci[j, k] <- 1
                if (res.ci[j, k] < 0)
                  res.ci[j, k] <- 0
              }
            }
          }
        }
        if (ntimes > 2) {
          if (auxi == 1) {
            for (j in 1:n.auxi) {
              be <- rep(0, n)
              be2 <- rep(0, n)
              for (i in 1:n) {
                be[i] <- Beran(ndata$Stime, 1 - ndata$event,
                               ndata[, covar], delta4, z.value[j], ndata$Stime[i],
                               kernel = window, bw = lbd2)
                be2[i] <- Beran(ndata[, 2 * ntimes - 3],
                                1 - ndata[, 2 * ntimes - 2], ndata[,
                                                                   covar], delta4, z.value[j], ndata[i,
                                                                                                     2 * ntimes - 3], kernel = window, bw = lbd2)
              }
              ifelse(method.weights == "NW", w1 <- NWW(ndata[,
                                                             covar], z.value[j], kernel = window, bw = lbd2),
                     w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value[j]))
              X <- data.frame(ndata[, 2 * (1:ntimes) -
                                      1])
              pos1 <- whichCS(X, x = x, lower.tail = lower.tail)
              pos2 <- which(ndata$Stime <= y & ndata$event ==
                              1)
              p1 <- intersect(pos1, pos2)
              pos3 <- which(ndata[, 2 * ntimes - 2] ==
                              1)
              p2 <- intersect(pos1, pos3)
              ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                              0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                       sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
              if (sum(w1[p2]) == 0)
                res.ci[j, k] <- NA
              if (is.na(res.ci[j, k]))
                count <- count + 1
              else {
                if (res.ci[j, k] > 1)
                  res.ci[j, k] <- 1
                if (res.ci[j, k] < 0)
                  res.ci[j, k] <- 0
              }
            }
          }
          if (auxi == 2) {

            be <- rep(0, n)
            be2 <- rep(0, n)
            for (i in 1:n) {
              be[i] <- Beran(ndata$Stime, 1 - ndata$event,
                             ndata[, covar], delta4, z.value, ndata$Stime[i],
                             kernel = window, bw = lbd2)
              be2[i] <- Beran(ndata[, 2 * ntimes - 3], 1 - ndata[, 2 * ntimes - 2], ndata[, covar],
                              delta4, z.value, ndata[i, 2 * ntimes - 3], kernel = window, bw = lbd2)
            }
            ifelse(method.weights == "NW", w1 <- NWW(ndata[, covar], z.value, kernel = window, bw = lbd2),
                   w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value))

            for (j in 1:n.auxi) {
              X <- data.frame(ndata[, 2 * (1:ntimes) - 1])
              pos1 <- whichCS(X, x = x, lower.tail = lower.tail)
              pos2 <- which(ndata$Stime <= y & ndata$event ==
                              1)
              p1 <- intersect(pos1, pos2)
              pos3 <- which(ndata[, 2 * ntimes - 2] ==
                              1)
              p2 <- intersect(pos1, pos3)
              ifelse(any(be[p1] == 0) | any(be2[p2] ==
                                              0), res.ci[j, k] <- NA, res.ci[j, k] <- 1 -
                       sum(w1[p1]/be[p1])/(sum(w1[p2]/be2[p2])))
              if (sum(w1[p2]) == 0)
                res.ci[j, k] <- NA
              if (is.na(res.ci[j, k]))
                count <- count + 1
              else {
                if (res.ci[j, k] > 1)
                  res.ci[j, k] <- 1
                if (res.ci[j, k] < 0)
                  res.ci[j, k] <- 0
              }
            }
          }
        }
        return(res.ci)
      }
      if (isTRUE(cluster)) {
        if (is.null(ncores)) {
          num_cores <- detectCores() - 1
        }
        else {
          num_cores <- ncores
        }
        registerDoParallel(cores = num_cores)
        on.exit(stopImplicitCluster())
        suppressMessages(res.ci <- foreach(i = 1:n.boot,
                                           .combine = cbind) %dorng% simplebootsurvIPCW(object,
                                                                                        n, n.auxi, auxi))
      }
      else {
        suppressMessages(res.ci <- foreach(i = 1:n.boot,
                                           .combine = cbind) %do% simplebootsurvIPCW(object,
                                                                                     n, n.auxi, auxi))
      }
      for (k in 1:n.auxi) {
        res.li[k] <- quantile(res.ci[k, ], (1 - conf.level)/2,
                              na.rm = TRUE)
        res.ls[k] <- quantile(res.ci[k, ], 1 - (1 - conf.level)/2,
                              na.rm = TRUE)
      }
      if (auxi == 1) {
        resu <- data.frame(cbind(y, res[!ii], res.li, res.ls))
        names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
      }
      if (auxi == 2) {
        resu <- data.frame(cbind(y, res[!ii], res.li, res.ls))
        names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
      }
    }

    text1 <- paste0("T", c(1:length(x)))
    text2 <- ifelse(lower.tail == TRUE, "<=", ">")
    text3 <- paste0(text1, text2, x, collapse = ",")
    callp <- paste("P(T>y|", text3, ",", z.name, "=", z.value, ")", sep = "")

    if (conf == FALSE & auxi == 1) {
      result <- list(est = resu, estimate = res[!ii], z.name = z.name,
                     z.value = z.value, x = x, y = y, bw = bw, window = window,
                     method.weights = method.weights, conf = conf, lbd = lbd2,
                     callp = callp)

    }

    if (conf == FALSE & auxi == 2) {
      result <- list(est = resu, estimate = res[!ii], z.name = z.name,
                     z.value = z.value, x = x, y = y, bw = bw, window = window,
                     method.weights = method.weights, conf = conf, lbd = lbd2,
                     callp = callp)
    }
    if (conf == TRUE & auxi == 1) {
      result <- list(est = resu, estimate = res[!ii], LCI = res.li,
                     UCI = res.ls, conf.level = conf.level, z.name = z.name,
                     z.value = z.value, bw = bw, window = window, method.weights = method.weights,
                     x = x, y = y, conf = conf, lbd = lbd2, callp = callp)

    }
    if (conf == TRUE & auxi == 2) {
      result <- list(est = resu, estimate = res[!ii], LCI = res.li,
                     UCI = res.ls, conf.level = conf.level, z.name = z.name,
                     z.value = z.value, bw = bw, window = window, method.weights = method.weights,
                     x = x, y = y, conf = conf, lbd = lbd2, callp = callp)

    }

    class(result) <- c("IPCW", "surv")
    return(invisible(result))
  }

