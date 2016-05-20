survCS <-
  function(time1, event1, Stime, event, ...)
  {
    if (missing(time1))
      stop("Argument 'time1' is missing, with no default")
    if (missing(event1))
      stop("Argument 'event1' is missing, with no default")
    if (missing(Stime))
      stop("Argument 'Stime' is missing, with no default")
    if (missing(event))
      stop("Argument 'event' is missing, with no default")

    data <- list(time1 = as.double(time1), event1 = as.integer(event1),
                 Stime = as.double(Stime), event = as.integer(event),
                 ...)
    datalen <- length(data)
    if (datalen > 4) {
      datanames <- names(data)
      for (i in 5:datalen) {
        if (!is.numeric(data[[i]]))
          stop("All additional arguments must be numeric")
        if (length(data[[i]]) != length(time1))
          stop("All additional arguments must have the same length as arguments 'time1', 'event1', 'Stime' and 'event'")
        if (datanames[i] == "")
          datanames[i] <- paste("covariate", i - 4, sep = ".")
        if (!is.double(data[[i]]))
          data[[i]] <- as.double(data[[i]])
      }
      names(data) <- datanames
    }
    attr(data, "row.names") <- as.integer(1:length(time1))
    attr(data, "class") <- "data.frame"
    #object <- vector(mode = "list", length = 1)
    #object <- na.omit(data)
    object <- list(data = na.omit(data))
    attr(object, "class") <- c("survCS", "surv")
    return(object)
  }
