#' Conditional survival probabilities based on the Kaplan-Meier weights,
#' Landmark approaches and Inverse probability of censoring weighted.
#'
#' Provides estimates for the conditional survival probabilities based on
#' Kaplan-Meier weighted estimators, the Landmark approaches and Inverse
#' probability of censoring weighted.
#'
#' Possible options for argument window are "gaussian", "epanechnikov",
#' "tricube", "boxcar", "triangular", "quartic" or "cosine".
#'
#' @param formula A formula object, which must have a \code{survCS} object as
#' the response on the left of the \code{~} operator and, if desired, a term on
#' the right. The term may be a qualitative or quantitative variable.  For a
#' single survival curve the right hand side should be ~ 1.
#' @param x Time or vector of times for the condional event(s).
#' @param y The total time for obtaining estimates for the conditional survival
#' probabilities.
#' @param lower.tail vector of logical values with the same size as 'x'. If 'x'
#' has dimension one and if \code{lower.tail = FALSE} (default), probabilities are \eqn{P(T
#' > y|T1 > x)} otherwise, \eqn{P(T > y|T1 <= x)}.  If the conditional event is
#' 2-dimensional, then, for example, given \eqn{x = c(x1, x2)} and \code{lower.tail =
#' c(TRUE, FALSE)} must be used to obtain probabilities \eqn{P(T > y|T1 <= x1, T2 >
#' x2)}. Multi-dimensional conditional events are introduced similarly.
#' @param method The method used to compute the conditional survival function.
#' Possible options are "LDM" and "KMW". Defaults to "LDM".
#' @param presmooth A logical value. If \code{TRUE}, the presmoothed landmark
#' estimator of the conditional survival function is computed. Only valid for
#' method = "LDM".
#' @param conf Provides pointwise confidence bands. Defaults to TRUE.
#' @param n.boot The number of bootstrap samples. Defaults to 200 samples.
#' @param data A data frame in which to interpret the variables named in the
#' \code{formula} argument.
#' @param conf.level Level of confidence. Defaults to 0.95 (corresponding to
#' 95\%).
#' @param z.value The value of the covariate on the right hand side of formula
#' at which the conditional survival probabilities are computed. For
#' quantitative covariates, i.e. of class integer and numeric.
#' @param bw A single numeric value to compute a kernel density bandwidth.  Use
#' "dpik" for the KernSmooth package based selector or "np" for the 'npudensbw'
#' function of the np package.
#' @param window A character string specifying the desired kernel.  See details
#' below for possible options. Defaults to "gaussian" where the gaussian
#' density kernel will be used.
#' @param method.weights A character string specifying the desired weights
#' method. Possible options are "NW" for the Nadaraya-Watson weights and "LL"
#' for local linear weights. Defaults to "NW".
#' @param cluster A logical value. If \code{TRUE} (default), the bootstrap
#' procedure for the confidence intervals is parallelized.  Note that there are
#' cases (e.g., a low number of bootstrap repetitions) that R will gain in
#' performance through serial computation. R takes time to distribute tasks
#' across the processors also it will need time for binding them all together
#' later on. Therefore, if the time for distributing and gathering pieces
#' together is greater than the time need for single-thread computing, it does
#' not worth parallelize.
#' @param ncores An integer value specifying the number of cores to be used in
#' the parallelized procedure. If \code{NULL} (default), the number of cores to
#' be used is equal to the number of cores of the machine - 1.
#' @param na.rm A logical value indicating whether NA values should be stripped
#' in the computation.
#' @return An object of class "survCS" and one of the following four classes:
#' "KMW", "LMD", "PLDM" and "IPCW". Objects are implemented as a list with
#' elements: \item{est}{data.frame with estimates of the conditional
#' probabilities.} \item{estimate}{Estimates of the conditional survival
#' probability.} \item{LCI}{The lower conditional survival probabilities of the
#' interval.} \item{UCI}{The upper conditional survival probabilities of the
#' interval.} \item{conf.level}{Level of confidence.} \item{y}{The total time
#' for obtaining the estimates of the conditional survival probabilities.}
#' \item{x}{The first time for obtaining the estimates of the conditional
#' survival probabilities.} \item{Nlevels}{The number of levels of the
#' covariate. Provides important information when the covariate at the right
#' hand side of formula is of class factor.} \item{conf}{logical; if FALSE
#' (default) the pointwise confidence bands are not given.} \item{callp}{The
#' expression of the estimated probability.} \item{levels}{The levels of the
#' qualitative covariate (if it is of class factor) on the right hand side of
#' formula.}
#'
#' @author Luis Meira-Machado and Marta Sestelo
#' @references L. Meira-Machado, M. Sestelo, and A. Goncalves (2016).
#' Nonparametric estimation of the survival function for ordered multivariate
#' failure time data: a comparative study. Biometrical Journal, 58(3),
#' 623--634.
#' @examples
#'
#'
#'    fit <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365, y = 730,
#'    data = colonCS, method = "KMW", conf = FALSE)
#'
#'    fit1 <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
#'    data = colonCS, method = "LDM", conf = FALSE)
#'
#'    fit2 <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
#'    data = colonCS, method = "LDM", lower.tail = TRUE, conf = FALSE)
#'
#'    fit3 <- survCOND(survCS(time1, event1, Stime, event) ~ 1, x = 365,
#'    y = c(730, 1095, 1460), data = colonCS, method = "LDM", presmooth = TRUE,
#'    lower.tail = TRUE, conf = TRUE, n.boot = 100, conf.level = 0.95,
#'    cluster = FALSE)
#'
#'    fit4 <- survCOND(survCS(time1, event1, Stime, event) ~ rx, x = 365,
#'    data = colonCS, method = "LDM", conf = FALSE)
#'
#'    fit5 <- survCOND(survCS(time1, event1, Stime, event) ~ factor(sex),
#'    x = 365, data = colonCS, method = "LDM", conf = FALSE)
#'
#'   \dontrun{
#'    fit6 <- survCOND(survCS(time1, event1, Stime, event) ~ age,
#'    x = 365, y = 730, z.value = 48, data = colonCS, conf = TRUE)
#'    }
#'
#'
#'
#'
#'
#' @export survCOND
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

