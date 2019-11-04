#Calulates Gini's mean difference.
#' @import parallelDist
#' @importfrom matrixStats rowDiffs
#' @importfrom RcppAlgos comboGeneral
#' @export
GMD <- function(Data) {
  len <- length(Data)
  cmb <- comboGeneral(Data,2, Parallel = TRUE)
  GMDn <- sum(abs(rowDiffs(cmb)))/(len^2*0.5)
  return(GMDn)
}