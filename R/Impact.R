# Calculates the Impact effect size measure.
#' @useDynLib(ImpactEffectsize, .registration = TRUE)
#' @importFrom caTools trapz
#' @importFrom methods hasArg
#' @importFrom matrixStats rowDiffs
#' @importFrom stats density var ks.test
#' @export
Impact <- function(Data, Cls, PlotIt = FALSE, pde = TRUE, col = c("red", "blue"),
                   meanLines = FALSE, medianLines = FALSE, ...) {
  # Convert Data and Cls to vectors
  Data <- as.vector(as.matrix(Data))
  Cls <- as.vector(as.matrix(Cls))
  
  # Check if Cls has exactly two distinct classes
  if (length(unique(Cls)) != 2) {
    stop("Impact: Cls has to contain exactly two distinct classes!")
  }
  
  # Check if Data and Cls have the same length
  if (length(Data) != length(Cls)) {
    stop("Impact: Data and Cls have different lengths!")
  }
  
  # Remove NAs
  if (length(which(is.na(Data))) > 0) {
    warning("NAs detected and removed.", call. = FALSE)
    Cls <- Cls[which(!is.na(Data))]
    Data <- Data[which(!is.na(Data))]
  }
  
  # Plot the densities if requested
  if (PlotIt == TRUE) {
    plot2Densities(Data = Data, Cls = Cls, col = col, pde = pde, meanLines = meanLines, medianLines = medianLines, ...)
  }
  
  # Initialize variables
  DirCT <- 1
  DirMorph <- 1
  UniqueCls <- sort(unique(Cls))
  ImpactX2X1 <- NA
  MorphDiff <- NA
  CTDiff <- NA
  GMDdata <- NA
  
  # Perform Kolmogorov-Smirnov test
  KStry <- try(suppressWarnings(stats::ks.test(Data[Cls == UniqueCls[1]], Data[Cls == UniqueCls[2]]))$p.value, TRUE)
  if (!inherits(KStry, "try-error")) {
    KSpval <- KStry
  } else {
    KSpval <- 1
  }
  
  # Calculate the Impact effect size
  if (stats::var(Data) == 0 || KSpval >= 0.05) {
    ImpactX2X1 <- 0
  } else {
    MedianCls1 <- c_median(Data[Cls == UniqueCls[1]])
    MedianCls2 <- c_median(Data[Cls == UniqueCls[2]])
    DeltaM <- MedianCls2 - MedianCls1
    DirCT <- ifelse(MedianCls2 < MedianCls1, -1, 1)
    GMDdata <- ClassGMD(Data, Cls)
    CTDiff <- abs(DeltaM) / GMDdata
    CTDiffWeight <- min(CTDiff, 2) / 2
    if ((stats::var(Data[Cls == UniqueCls[1]]) == 0 ||
         stats::var(Data[Cls == UniqueCls[2]]) == 0) && stats::var(Data) > 0) {
      MorphDiff <- 0
    } else {
      set.seed(42)
      PDEKernelsTry <- try(ParetoDensityEstimationIE(Data), TRUE)
      if (!inherits(PDEKernelsTry, "try-error")) {
        PDEKernels <- PDEKernelsTry$kernels
        pdeX1Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[1]], kernels = PDEKernels), TRUE)
        pdeX2Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[2]], kernels = PDEKernels), TRUE)
      }
      MorphDiff <- 0
      if (!inherits(PDEKernelsTry, "try-error") &&
          !inherits(pdeX1Try, "try-error") && !inherits(pdeX2Try, "try-error")) {
        pdeX1 <- pdeX1Try$paretoDensity
        pdeX2 <- pdeX2Try$paretoDensity
        pdeDiff <- abs(matrixStats::rowDiffs(cbind(pdeX2, pdeX1)))
        Momentum1 <- (sum(sign(Data[Cls == UniqueCls[1]]) * log10(abs(Data[Cls == UniqueCls[1]]) + 1))) /
          (sign(length(Data[Cls == UniqueCls[1]])) * log10(abs(length(Data[Cls == UniqueCls[1]]) + 1)))
        Momentum2 <- (sum(sign(Data[Cls == UniqueCls[2]]) * log10(abs(Data[Cls == UniqueCls[2]]) + 1))) /
          (sign(length(Data[Cls == UniqueCls[2]])) * log10(abs(length(Data[Cls == UniqueCls[2]]) + 1)))
        DirMorph <- ifelse(Momentum2 < Momentum1, -1, 1)
        if (length(PDEKernels) == length(pdeDiff)) {
          MorphDiff <- caTools::trapz(PDEKernels, pdeDiff)
        }
      }
    }
    ImpactX2X1 <- CTDiffWeight * (DirCT * abs(CTDiff)) + (1 - CTDiffWeight) * (DirMorph * abs(MorphDiff))
  }
  
  return(list("Impact" = ImpactX2X1, "MorphDiff" = MorphDiff, "CTDiff" = CTDiff))
}
