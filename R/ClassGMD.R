# Pools two values of Gini's mean difference calculated for separate groups.
#' @importFrom stats var
ClassGMD <- function(Data, Cls) {
  # Set a maximum number of data points
  maxPoints <- 10000
  
  # If the data has more than the maximum points, sample the data
  if (length(Data) > maxPoints) {
    df <- data.frame(cbind(Data, Cls))
    table(df$Cls)
    dfsplit <- split(df, list(df$Cls))
    set.seed(42)
    samples <- lapply(dfsplit, function(x) {
      x[sample(1:nrow(x), maxPoints, FALSE), ]
    })
    out <- do.call(rbind, samples)
    table(out$Cls)
    Data <- as.vector(out$Data)
    Cls <- as.vector(out$Cls)
  }
  
  # Calculate the Gini's mean difference (GMD)
  if (stats::var(Data) == 0) {
    GMDn <- 1e-7
  } else if ((stats::var(Data[Cls == sort(unique(Cls))[1]]) == 0 |
              stats::var(Data[Cls == sort(unique(Cls))[2]]) == 0) &
             stats::var(Data) > 0) {
    GMDn <- c_gmd(Data)
  } else {
    GMD1 <- c_gmd(Data[Cls == sort(unique(Cls))[1]])
    GMD2 <- c_gmd(Data[Cls == sort(unique(Cls))[2]])
    GMDn <- sqrt((GMD1^2 + GMD2^2) / 2)
  }
  
  return(GMDn)
}
