#Calculates the Impact efefct size measure.
#' @importFrom caTools trapz
#' @importFrom AdaptGauss ParetoDensityEstimation
#' @importFrom methods hasArg
#' @importFrom stats density median var
#' @export
Impact <- function(Data, Cls, PlotIt = FALSE, pde = TRUE, 
  col = c("red","blue"), meanLines = FALSE, medianLines = FALSE, ...){
  if(length(unique(Cls))!=2) stop("Impact: Cls has to contain exactly two distinct classes!")
  if(length(Data) != length(Cls)) stop("Impact: Data and Cls have different lengths!")
  if(hasArg(PlotIt) == TRUE & PlotIt == TRUE) plot2Densities(Data,Cls, ...)
  ImpactX2X1 <- NA; MorphDiff = NA; CTDiff = NA; GMDdata = NA 
  if(var(Data)==0) ImpactX2X1 = 0 else {
    DeltaM = median(Data[Cls==unique(Cls)[2]]) -  median(Data[Cls == unique(Cls)[1]])
    Direction <- ifelse(mean(Data[Cls==unique(Cls)[2]]) < mean(Data[Cls == unique(Cls)[1]]),-1,1)          
    GMDdata <- ClassGMD(Data,Cls)
    CTDiff <- DeltaM/GMDdata
    if((var(Data[Cls==unique(Cls)[1]]) == 0 | var(Data[Cls==unique(Cls)[2]]) == 0) & var(Data) > 0) {
      Overlap <- 0
      MorphDiff <- 0
    } else {
      PDEKernels <- AdaptGauss::ParetoDensityEstimation(Data)$kernels
      pdeX1 <- ParetoDensityEstimation(Data = Data[Cls==unique(Cls)[1]], kernels = PDEKernels)$paretoDensity
      pdeX2 <- ParetoDensityEstimation(Data = Data[Cls==unique(Cls)[2]], kernels = PDEKernels)$paretoDensity
      pdeDiff <- abs(rowDiffs(cbind(pdeX2,pdeX1)))
      MorphDiff <- Direction*trapz(PDEKernels,pdeDiff)
    }
    ImpactX2X1 <- CTDiff + MorphDiff
  }
  return(list("Impact" = ImpactX2X1, "MorphDiff" = MorphDiff, "CTDiff" = CTDiff, "GMDdata" = GMDdata))
}
