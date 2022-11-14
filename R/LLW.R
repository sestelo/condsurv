#' Local linear weights.
#' 
#' Computes local linear weights based on Kernel smoothing.
#' 
#' Possible options for argument window are "gaussian", "epanechnikov",
#' "tricube", "boxcar", "triangular", "quartic" or "cosine".
#' 
#' @param x Covariate values for obtaining estimates for the conditional
#' probabilities. If missing, unconditioned probabilities will be computed.
#' @param kernel A character string specifying the desired kernel. See details
#' below for possible options. Defaults to "gaussian" where the gaussian
#' density kernel will be used.
#' @param bw A single numeric value to compute a kernel density bandwidth.
#' @param t1 Covariate value to compute the weight at.
#' @return A vector with local linear weights.
#' @author Luis Meira-Machado and Marta Sestelo
#' @seealso \code{\link{NWW}}
#' @examples
#' 
#' LLW(x = colonCS$age, bw = 3, t1 = 60)
#' 
LLW <-
  function(x, kernel = "gaussian", bw, t1){
    nobs <- length(x)
    return(.C("LLWeightsKernel",as.double(x), as.integer(nobs),
              as.double(t1), as.double(bw), as.character(kernel),
              w = double(nobs) )$w)}
