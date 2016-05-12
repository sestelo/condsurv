plot.surv <-
function(x = object, y = NULL, type = NULL, ...) {
  
  object <- x
  
  if (class(object)[1] != "KMW" & class(object)[1] != "LDM" & 
        class(object)[1] != "PLDM" & 
        class(object)[1] != "IPCW"){
    
    stop("The argumment 'Object' must be of one of the following classes
         'KMW', 'LDM', 'PLDM', 'IPCW'")
  }
  
  ob <- object$est
  
  if(is.null(type)){
    type <- "l"
  }
  
  if (object$conf == FALSE) { plot(ob[, 1], ob[, 2], type = type, ...) }
  
  if (object$conf == TRUE) {
    matplot(x = ob[, 1], y = ob[, 2:4], type = type, ...)
    #oask <- devAskNewPage(TRUE)
    #on.exit(devAskNewPage(oask))
  }

}
