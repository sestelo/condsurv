#' Nadaraya-Watson weights.
#' 
#' Computes the Nadaraya-Watson weights.
#' 
#' Possible options for argument window are "gaussian", "epanechnikov",
#' "tricube", "boxcar", "triangular", "quartic" or "cosine".
#' 
#' @param covariate Covariate values for obtaining weights.
#' @param x Covariate value to compute the weight at.
#' @param kernel A character string specifying the desired kernel. See details
#' below for possible options. Defaults to "gaussian" where the gaussian
#' density kernel will be used.
#' @param bw A single numeric value to compute a kernel density bandwidth.
#' @return A vector with Nadaraya-Watson weights.
#' @author Luis Meira-Machado and Marta Sestelo
#' @seealso \code{\link{LLW}}
#' @examples
#' 
#' NWW(covariate = colonCS$age, x=40, kernel = "gaussian", bw = 3)
#' 
NWW <-
  function(covariate, x, kernel="gaussian", bw) {
    spa <- NULL
    len <- length(covariate)
    listg <- .C( "NWWeightsKernel", as.double(covariate), as.integer(len),
                 as.double(x), as.double(bw), as.character(kernel),
                 weight = double(len) )
    return(listg$weight)}
