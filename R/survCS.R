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

    data <- data.frame(as.double(time1), as.integer(event1),
                       as.double(Stime), as.integer(event),
                       ...)

    names(data)[1:4] <- c(substitute(time1), substitute(event1),
                          substitute(Stime), substitute(event))

    data <- na.omit(data)

    datalen <- length(data)
    odd <- 1 - datalen %% 2
    ntimes <- datalen %/% 2

    for (i in 1:ntimes) {
      if (i != ntimes) {if(any(data[,2*i+1] < data[,2*i-1]))
        stop(paste("Values of the time variable", names(data)[2*i+1],
                   "can not be less than those for", names(data)[2*i-1])) }
      if (any(data[,2*i] != 0 & data[,2*i] != 1))
        stop("The status indicator variable must be 0 if the total time is censored and 1 otherwise.")
    }

    datanames <- names(data)

    for (i in 1:ntimes) {
      if (!is.numeric(data[,2*i-1]))
        stop("All arguments must be numeric")
      if (length(data[,2*i-1]) != length(data[,1]))
        stop("All arguments must have the same length")
      if (length(data[,2*i]) != length(data[,1]))
        stop("All arguments must have the same length")
      datanames[2*i-1] <- paste("time", i, sep = "")
      datanames[2*i] <- paste("event", i, sep = "")
    }
    datanames[2*ntimes-1] <- "Stime"
    datanames[2*ntimes] <- "event"
    names(data) <- datanames

    if (odd == 0) stop("Incorrect number of variables in 'survCS' ")


    attr(data, "row.names") <- as.integer(1:length(time1))
    attr(data, "class") <- "data.frame"
    object <- list(data = na.omit(data))
    attr(object, "class") <- c("survCS", "surv")
    return(object)
  }
