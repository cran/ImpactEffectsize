#Calculates the Impact efefct size measure.
#' @importFrom pracma trapz
#' @importFrom DataVisualizations ParetoDensityEstimation
#' @importFrom methods hasArg
#' @importFrom matrixStats rowDiffs
#' @importFrom stats density median var ks.test
#' @export
Impact <- function(Data, Cls, PlotIt = FALSE, pde = TRUE, 
                   col = c("red","blue"), meanLines = FALSE, medianLines = FALSE, ...){
  Data = as.vector(as.matrix(Data))
  Cls = as.vector(as.matrix(Cls))
  if(length(sort(unique(Cls)))!=2) stop("Impact: Cls has to contain exactly two distinct classes!")
  if(length(Data) != length(Cls)) stop("Impact: Data and Cls have different lengths!")
  if(hasArg(PlotIt) == TRUE & PlotIt == TRUE) plot2Densities(Data,Cls, ...)
  if(length(which(is.na(Data))) > 0) {
    warning("NAs detected and removed.", call. = FALSE)
    Cls <- Cls[which(!is.na(Data))]
    Data <- Data[which(!is.na(Data))]
  }
  DirCT<- 1
  DirMorph<- 1
  ImpactX2X1 <- NA; MorphDiff = NA; CTDiff = NA; GMDdata = NA 
  if(var(Data)==0 | ks.test(Data[Cls==sort(unique(Cls))[1]], Data[Cls==sort(unique(Cls))[2]])$p.value >= 0.05) ImpactX2X1 = 0 else {
    DeltaM <- median(Data[Cls==sort(unique(Cls))[2]]) -  median(Data[Cls == sort(unique(Cls))[1]])
    DirCT<- ifelse(median(Data[Cls==sort(unique(Cls))[2]]) < median(Data[Cls == sort(unique(Cls))[1]]),-1,1)          
    GMDdata <- ClassGMD(Data,Cls)
    CTDiff <- abs(DeltaM)/GMDdata
    CTDiffWeight <- min(CTDiff,2)/2
    if((var(Data[Cls==sort(unique(Cls))[1]]) == 0 | var(Data[Cls==sort(unique(Cls))[2]]) == 0) & var(Data) > 0) {
      Overlap <- 0
      MorphDiff <- 0
    } else {
      errorF0=try(DataVisualizations::ParetoDensityEstimation(Data),TRUE)
      if(class(errorF0) != "try-error") {
        set.seed(42)
        PDEKernels <- DataVisualizations::ParetoDensityEstimation(Data)$kernels
        errorF1=try(ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[1]], kernels = PDEKernels),TRUE)
        errorF2=try(ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[2]], kernels = PDEKernels),TRUE)
      }
      MorphDiff <- 0
      if(class(errorF0) != "try-error") {
        if (class(errorF1) != "try-error" & class(errorF2) != "try-error") {
          set.seed(42)
          pdeX1 <- ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[1]], kernels = PDEKernels)$paretoDensity
          pdeX2 <- ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[2]], kernels = PDEKernels)$paretoDensity
          pdeDiff <- abs(rowDiffs(cbind(pdeX2,pdeX1)))
          Momentum1 = (sum(sign(Data[Cls==sort(unique(Cls))[1]])*log10(abs(Data[Cls==sort(unique(Cls))[1]])+1)))/(sign(length(Data[Cls==sort(unique(Cls))[1]]))*log10(abs(length(Data[Cls==sort(unique(Cls))[1]])+1)))
          Momentum2 = (sum(sign(Data[Cls==sort(unique(Cls))[2]])*log10(abs(Data[Cls==sort(unique(Cls))[2]])+1)))/(sign(length(Data[Cls==sort(unique(Cls))[2]]))*log10(abs(length(Data[Cls==sort(unique(Cls))[2]])+1)))
          DirMorph<- ifelse(Momentum2 < Momentum1,-1,1)
          if (length(PDEKernels) == length(pdeDiff)) MorphDiff <- trapz(PDEKernels,pdeDiff)
        }
      }
    }
    ImpactX2X1 <- CTDiffWeight * (DirCT* abs(CTDiff)) + (1 - CTDiffWeight) * (DirMorph* abs(MorphDiff))
  }
  return(list("Impact" = ImpactX2X1, "MorphDiff" = MorphDiff, "CTDiff" = CTDiff, "GMDdata" = GMDdata))
}
