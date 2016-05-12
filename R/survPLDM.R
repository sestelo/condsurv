survPLDM <-
function(object, x, y, conf = FALSE, n.boot = 1000, conf.level = 0.95, lower.tail = FALSE)
{
    if (missing(object))
        stop("Argument 'object' is missing, with no default")
    if (missing(x))
        x <- 0
    if (missing(y))
        y <- max(object[[1]]$Stime)
    if (length(x) > 1) stop("Length of 'x' must be 1")

    y <- y[ y>=x ]
    y <- sort(unique(y))
    p1 <- which(object[[1]]$time1 > x)
    if(lower.tail == TRUE) p1 <- which(object[[1]]$time1 <= x)
    G0 <- PKMW(object[[1]]$Stime[p1], object[[1]]$event[p1])
    t_2 <- object[[1]]$Stime[p1]
    res <- rep(0, length(y))
    for (k in 1: length(y)) {
                p2 <- which(t_2 <= y[k])
                res[k] <- 1 - sum(G0[p2])
                            }
    res.li <- rep(0,length(y))
    res.ls <- rep(0,length(y))
    resu <- data.frame(cbind(y, res))
    names(resu) <- c("y", "estimate")


    if (conf == TRUE)
           {
           res.ci <- matrix(0, nrow = length(y), ncol = n.boot)
           for (j in 1:n.boot)
                   {
                   n <- dim(object[[1]])[1]
                   xx <- sample.int(n, size = n, replace = TRUE)
                   ndata <- object[[1]][xx,]
                   p1 <- which(ndata$time1 > x)
   if (lower.tail == TRUE) p1 <- which(ndata$time1 <= x)
                   G0 <- PKMW(ndata$Stime[p1], ndata$event[p1])
                   t_2 <- ndata$Stime[p1]
                   for (k in 1: length(y)) {
                                  p2 <- which(t_2 <= y[k])
                                  res.ci[k,j] <- 1 - sum(G0[p2])
                                                 }
                   }
                   for (k in 1: length(y)) {
                                  res.li[k] <- quantile(res.ci[k,], ( 1 - conf.level) / 2)
                                  res.ls[k] <- quantile(res.ci[k,], 1 - (1 - conf.level) / 2)
                                           }
           if (length(y)==1 & lower.tail == FALSE) cat("S(T>",y,"|T1>",x,") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
   if (length(y)==1 & lower.tail == TRUE) cat("S(T>",y,"|T1<=",x,") = ", res,"  ", conf.level*100,"%CI: ", res.li, "-", res.ls, sep="", "\n")
           if (length(y)>1) {
                         resu <- data.frame(cbind(resu, res.li, res.ls))
                         names(resu) <- c("y", "estimate", "LCI", "UCI")
                         if (lower.tail == FALSE) cat("Estimates of S(T>y,|T1>",x,")",sep="","\n")
 if (lower.tail == TRUE) cat("Estimates of S(T>y,|T1<=",x,")",sep="","\n")
                         print(resu)
                            }
           }

    if(conf == FALSE) {
                  result <- list(est = resu, estimate = res, y = y, x = x, conf = conf)
                  if (length(y) == 1 & lower.tail == FALSE) cat("S(T>",y,"|T1>",x,") = ", res, sep="", "\n")
  if (length(y) == 1 & lower.tail == TRUE) cat("S(T>",y,"|T1<=",x,") = ", res, sep="", "\n")
                  if (length(y)>1) {
                         if (lower.tail  == FALSE) cat("Estimates of S(T>y,|T1>",x,")",sep="","\n")
 if (lower.tail  == TRUE) cat("Estimates of S(T>y,|T1<=",x,")",sep="","\n")
                         print(resu)
                                   }
                    }
    if(conf==TRUE) {
                  result <- list(est = resu, estimate = res, LCI = res.li, UCI = res.ls, conf.level = conf.level, y = y, x = x, conf = conf)
                   }
    class(result) <- c("PLDM", "surv")
    return(invisible(result))
}
