#Pools two values of Gini's mean difference calulated for separate groups.
#' @importFrom stats var
ClassGMD <- function(Data, Cls) {
  if(var(Data) == 0) GMDn <- 1e-7 
  else if ((var(Data[Cls==unique(Cls)[1]]) == 0 | var(Data[Cls==unique(Cls)[2]]) == 0) & var(Data) > 0) GMDn <- GMD(Data) else {
    GMD1 <- GMD(Data[Cls==unique(Cls)[1]])
    GMD2 <- GMD(Data[Cls==unique(Cls)[2]])
    GMDn <- sqrt(GMD1^2+GMD2^2)
  } 
  return(GMDn)
}
