whichCS <- function(times, x, lower.tail)
{
  A <- matrix(NA, ncol = length(x), nrow = dim(times)[1])
  B <- times
  for (j in 1:length(x)) {
    ifelse(lower.tail[j] == TRUE, A[,j] <- B[,j] <= x[j], A[,j] <- B[,j] > x[j])
  }
  pos <- which(rowSums(A) == length(x))
  return(pos)
}

