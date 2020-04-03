#Pools two values of Gini's mean difference calculated for separate groups.
#' @importFrom stats var
#' @importFrom Hmisc GiniMd
ClassGMD <- function(Data, Cls) {
  if(var(Data) == 0) GMDn <- 1e-7 
  else if ((var(Data[Cls==sort(unique(Cls))[1]]) == 0 | var(Data[Cls==sort(unique(Cls))[2]]) == 0) & var(Data) > 0) GMDn <- GiniMd(Data) else {
    GMD1 <- GiniMd(Data[Cls==sort(unique(Cls))[1]])
    GMD2 <- GiniMd(Data[Cls==sort(unique(Cls))[2]])
    GMDn <- sqrt((GMD1^2+GMD2^2)/2)
  } 
  return(GMDn)
}
