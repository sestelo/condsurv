print.surv <- function(x, ...){

  if (inherits(x, "surv")) {
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:\n")
  print(class(x)[1])

  }else{
    stop("Argument x must be either surv object.")
}
}
