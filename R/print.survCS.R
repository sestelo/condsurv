print.survCS <- function(x, ...){

  if (inherits(x, "survCS")) {

    if(class(x)[1] == "data.frame") {

      print(x)

    }else{


      cat("Call:\n")
      print(x$call)
      cat("\nMethod:\n")

      if(class(x)[1] == "KMW") method <- "Kaplan-Meier weights"
      if(class(x)[1] == "LDM") method <- "Landmark approach"
      if(class(x)[1] == "PLDM") method <- "Presmoothed Landmark approach"
      if(class(x)[1] == "IPCW") method <- "Inverse Probability of Censoring Weighting"

      print(method)
    }

  }else{
    stop("Argument x must be either survCS object.")
  }
}
