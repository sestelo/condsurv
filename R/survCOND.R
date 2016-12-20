survCOND <- function(formula, x, y, lower.tail = FALSE, method = "LDM",
                     presmooth = FALSE, conf = TRUE, n.boot = 200, data,
                     conf.level = 0.95, z.value, bw = "dpik", window = "gaussian",
                     method.weights = "NW", cluster = FALSE, ncores = NULL, na.rm = TRUE)
{

  if (missing(formula)) stop("A formula argument is required")
  if (missing(x)) stop("argument 'x' is missing, with no default")
  if (length(x) != length(lower.tail)) stop("Arguments 'x' and 'lower.tail' must have the same length")
  if (any(x < 0)) stop("'x' values should be nonnegative")
  if(is.unsorted(x)) stop("'x' values should be sorted")

  fmla <- eval(formula, parent.frame())
  Terms <- terms(fmla)
  mf <- Terms[[2]]

  object <- with(data = data, eval(mf))
  if (!inherits(object, "survCS")) stop("Response must be a survCS object")

  object <- list(data = object) # new since survCS doesn't return a list

  obj_data <- object[[1]]

  X <- Terms[[3]]
  if(length(attr(terms(formula),"term.labels")) > 1) stop("only one covariate is supported")
  Class <- class(with(data=data, eval(Terms[[3]])))
  if (Class != "numeric" & Class != "integer" & Class != "factor") stop("covariate must be one of the following classes 'numeric', 'integer' or 'factor'")
  xval <- with(data=data, eval(X))
  lencov <- length(xval)
  lencov2 <- length(attr(terms(formula),"term.labels"))
  if(lencov2 != 0) Xval <- with(data=data, eval(Terms[[3]]))
  if(lencov != dim(obj_data)[1] & lencov2 != 0) stop("length of the covariate does not match")

  lenc <- dim(object[[1]])[2]
  ntimes <- lenc%/%2
  if (length(x) != ntimes-1) {
    cat("The number of consecutive event times in 'survCS' is", ntimes, ". The length of 'x' should be", ntimes-1,"\n")
    stop("The length of 'x' is not supported for the selected 'object'")
  }


  if(missing(y)) {
    y <- sort(unique(obj_data[,dim(obj_data)[2]-1]))
    if(any(!lower.tail)) {
      pos <- max(which(!lower.tail))
      y <- y[y>=x[pos]]
    }
  }
  else {
    if(any(!lower.tail)){
      pos <- max(which(!lower.tail))
      if(any(y < x[pos])) stop("'y' values should be equal or greater than all values in 'x'")
    }
  }

  if (any(y < 0)) stop("'y' values should be nonnegative")
  y <- sort(unique(y))

  text1 <- paste0("T", c(1:length(x)))
  text2 <- ifelse(lower.tail == TRUE, "<=", ">")
  text3 <- paste0(text1, text2, x, collapse = ",")
  callp <- paste("P(T>y|", text3, ")", sep = "")


  if (!(method %in% c("LDM", "KMW"))){
    stop("Possible methods are 'LDM' and 'KMW'." )
  }

  if (method == "KMW" & presmooth == TRUE)  {
    warning("Argument 'presmooth' is not used by the 'KMW' method.")
  }

  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW    #remover argumento z.name?

    obj1 <- object
    obj1[[1]] <- cbind(obj1[[1]], xval)
    obj1[[1]] <- na.omit(obj1[[1]])
    colnames(obj1[[1]]) <- c(colnames(object[[1]]), attr(terms(formula),"term.labels"))

    res <- survIPCW(object = obj1, x = x, y =y, z.name = attr(terms(formula),"term.labels"),
                    z.value = z.value, bw = bw, window = window,
                    method.weights = method.weights, conf = conf, n.boot = n.boot,
                    conf.level = conf.level, lower.tail = lower.tail, cluster = cluster, ncores = ncores)

    class(res) <- c("IPCW", "survCS")

    if (length(y) == 1) resu <- res$est
    if (length(y) > 1 & conf == TRUE) {
      resu <- data.frame(cbind(res$y, res$estimate, res$LCI, res$UCI))
      names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
    }

    if (length(y)>1 & conf == FALSE) {
      resu <- data.frame(cbind(res$y, res$estimate))
      names(resu) <- c("y","estimate")
    }

    if(conf == TRUE) result <- list(est=resu, estimates=res$estimate, LCI=res$LCI, UCI=res$UCI, conf.level=conf.level, y=y, x=x, Nlevels=1, conf=conf, callp = res$callp,
                                    levels = NULL, call = match.call())
    if(conf == FALSE) result <- list(est=resu, estimates=res$estimate, conf.level=conf.level, y=y, x=x, Nlevels=1, conf=conf, callp = res$callp, levels = NULL, call = match.call())
    class(result) <- class(res)
    return(invisible(result))
  }



  if (length(attr(terms(formula), "term.labels")) == 0) {  #LDM/PLMD/KMW

    if (method == "LDM" & presmooth == FALSE)  {
      res <- survLDM(object = object, x = x, y = y, conf = conf,
                     n.boot = n.boot, conf.level = conf.level,
                     lower.tail = lower.tail, cluster = cluster,
                     ncores = ncores)
      class(res) <- c("LDM", "survCS")
    }

    if (method == "LDM" & presmooth == TRUE) {
      res <- suppressWarnings(survPLDM(object = object, x = x, y = y,
                                       conf = conf, n.boot = n.boot,
                                       conf.level = conf.level,
                                       lower.tail = lower.tail,
                                       cluster = cluster, ncores = ncores))
      class(res) <- c("PLDM", "survCS")
    }

    if (method == "KMW")  {
      res <- survKMW(object = object, x = x, y = y, conf = conf,
                     n.boot = n.boot, conf.level = conf.level,
                     lower.tail = lower.tail, cluster = cluster,
                     ncores = ncores, na.rm = na.rm)
      class(res) <- c("KMW", "survCS")
    }


    #if (length(y) == 1) resu <- res$est
    if (length(y) >= 1 & conf == TRUE) {
      resu <- data.frame(cbind(res$y, res$estimate, res$LCI, res$UCI))
      names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
    }

    if (length(y)>=1 & conf == FALSE) {
      resu <- data.frame(cbind(res$y, res$estimate))
      names(resu) <- c("y","estimate")
    }



    if(conf == TRUE) result <- list(est=resu, estimates=res$estimate, LCI=res$LCI, UCI=res$UCI, conf.level=conf.level, y=y, x=x, Nlevels=1, conf=conf, callp = callp,
                                    levels = NULL, call = match.call())
    if(conf == FALSE) result <- list(est=resu, estimates=res$estimate, conf.level=conf.level, y=y, x=x, Nlevels=1, conf=conf, callp = callp, levels = NULL, call = match.call())
    class(result) <- class(res)
    return(invisible(result))
  }

  if (length(attr(terms(formula),"term.labels")) > 0 & Class == "factor") {  #LDM/PLMD/KMW by levels of the covariate
    x.nlevels <- nlevels(with(data=data, eval(formula[[3]])))
    levels <- levels(with(data=data, eval(formula[[3]])))

    estim <- c()
    LCI <- c()
    UCI <- c()

    for (k in 1:x.nlevels) {
      v.level <- levels(with(data=data, eval(formula[[3]])))[k]
      p<- which(Xval == v.level)
      obj<- object
      obj$data <- object$data[p,]

      if (method == "LDM" & presmooth == FALSE)  {
        res <- survLDM(object = obj, x = x, y = y, conf = conf,
                       n.boot = n.boot, conf.level = conf.level,
                       lower.tail = lower.tail, cluster = cluster,
                       ncores = ncores)
        class(res) <- c("LDM", "survCS")
      }

      if (method == "LDM" & presmooth == TRUE) {
        res <- suppressWarnings(survPLDM(object = obj, x = x, y = y,
                                         conf = conf, n.boot = n.boot,
                                         conf.level = conf.level,
                                         lower.tail = lower.tail,
                                         cluster = cluster, ncores = ncores))
        class(res) <- c("PLDM", "survCS")
      }

      if (method == "KMW")  {
        res <- survKMW(object = obj, x = x, y = y, conf = conf,
                       n.boot = n.boot, conf.level = conf.level,
                       lower.tail = lower.tail, cluster = cluster,
                       ncores = ncores, na.rm = na.rm)
        class(res) <- c("KMW", "survCS")
      }

      if (conf == TRUE) {
        resu <- data.frame(cbind(res$y, res$estimate, res$LCI, res$UCI))
        names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
      }

      if (conf == FALSE) {
        resu <- data.frame(cbind(res$y, res$estimate))
        names(resu) <- c("y","estimate")
      }

      estim[[paste(levels[k])]] <- resu
    }



    if (conf == TRUE){
      result <- list(est=estim, conf.level=conf.level, y=y, x=x, Nlevels=x.nlevels, conf=conf, callp = callp, formula = formula, levels = levels, call = match.call())
      class(result) <- class(res)
      return(invisible(result))
    }

    if (conf == FALSE){
      result <- list(est=estim, y=y, x=x, Nlevels=x.nlevels, conf=conf, callp = callp, formula = formula, levels = levels, call = match.call())
      class(result) <- class(res)
      return(invisible(result))
    }
  }
}

