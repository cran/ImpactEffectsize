#' Calculate the Impact effect size measure between two groups.
#'
#' This function computes an effect size measure based on central tendency and morphological differences between two groups.
#' It supports Pareto Density Estimation (PDE) for morphology and falls back to standard kernel density if PDE fails.
#' Optionally, it plots the density estimates.
#'
#' @param Data Numeric vector of data points.
#' @param Cls Grouping vector or factor with exactly two distinct classes.
#' @param PlotIt Logical; if TRUE, plots the density estimates.
#' @param pde Logical; if TRUE, attempts Pareto density estimation, else uses kernel density.
#' @param col Character vector of length 2; colors for the two groups in plots.
#' @param meanLines Logical; if TRUE, adds mean lines to the plot.
#' @param medianLines Logical; if TRUE, adds median lines to the plot.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return A list containing:
#' \describe{
#'   \item{Impact}{Numeric; the combined impact effect size measure.}
#'   \item{MorphDiff}{Numeric; morphological difference component.}
#'   \item{CTDiff}{Numeric; central tendency difference component.}
#'   \item{density_df}{Data frame; columns PDEKernels, pde_Cls1, pde_Cls2 containing density data used for plotting.}
#' }
#' 
#' @importFrom caTools trapz
#' @importFrom methods hasArg
#' @importFrom matrixStats rowDiffs
#' @importFrom stats density var ks.test
#' @importFrom withr with_seed
#' @export
Impact <- function(Data, Cls, PlotIt = FALSE, pde = TRUE, col = c("red", "blue"),
                   meanLines = FALSE, medianLines = FALSE, ...) {
  
  # Convert Data and Cls to vectors
  Data <- as.vector(as.matrix(Data))
  Cls <- as.vector(as.matrix(Cls))
  
  # Validate Cls
  if (length(unique(Cls)) != 2) {
    stop("Impact: Cls has to contain exactly two distinct classes!")
  }
  
  # Validate lengths
  if (length(Data) != length(Cls)) {
    stop("Impact: Data and Cls have different lengths!")
  }
  
  # Remove NAs from Data and corresponding Cls entries
  if (anyNA(Data)) {
    warning("NAs detected and removed.", call. = FALSE)
    keep <- !is.na(Data)
    Data <- Data[keep]
    Cls  <- Cls[keep]
  }
  
  # Initialize variables
  DirCT <- 1
  DirMorph <- 1
  UniqueCls <- sort(unique(Cls))
  ImpactX2X1 <- NA
  MorphDiff <- 0
  CTDiff <- NA
  GMDdata <- NA
  PDEKernels <- NULL
  pdeX1 <- NULL
  pdeX2 <- NULL
  
  # Perform Kolmogorov-Smirnov test for equality of distributions
  KStry <- try(suppressWarnings(stats::ks.test(Data[Cls == UniqueCls[1]], Data[Cls == UniqueCls[2]]))$p.value, TRUE)
  KSpval <- if (!inherits(KStry, "try-error")) KStry else 1
  
  # If no variability or too few points, Impact set to zero and no PDE possible
  if (length(unique(Data)) <= 1 || length(Data) < 2) {
    ImpactX2X1 <- 0
    if (PlotIt) {
      plot2Densities(Data, Cls, pde = FALSE, meanLines = meanLines, medianLines = medianLines, col = col, ...)
    }
  } else {
    # Calculate medians and delta
    MedianCls1 <- c_median(Data[Cls == UniqueCls[1]])
    MedianCls2 <- c_median(Data[Cls == UniqueCls[2]])
    DeltaM <- MedianCls2 - MedianCls1
    DirCT <- ifelse(MedianCls2 < MedianCls1, -1, 1)
    
    # Generalized Mean Difference (GMD)
    GMDdata <- ClassGMD(Data, Cls)
    CTDiff <- abs(DeltaM) / GMDdata
    CTDiffWeight <- min(CTDiff, 2) / 2
    
    # Check zero variance in groups, fallback if so
    if ((stats::var(Data[Cls == UniqueCls[1]], na.rm = TRUE) == 0 ||
         stats::var(Data[Cls == UniqueCls[2]], na.rm = TRUE) == 0) &&
        stats::var(Data, na.rm = TRUE) > 0) {
      MorphDiff <- 0
      if (PlotIt) {
        plot2Densities(Data, Cls, pde = FALSE, meanLines = meanLines, medianLines = medianLines, col = col, ...)
      }
    } else {
      # Try Pareto Density Estimation with fixed seed for reproducibility
      withr::with_seed(42, {
        PDEKernelsTry <- try(ParetoDensityEstimationIE(Data), TRUE)
      })
      
      # Attempt group-wise PDE if global PDE succeeds
      if (!inherits(PDEKernelsTry, "try-error")) {
        PDEKernels <- PDEKernelsTry$kernels
        pdeX1Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[1]], kernels = PDEKernels), TRUE)
        pdeX2Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[2]], kernels = PDEKernels), TRUE)
      } else {
        pdeX1Try <- pdeX2Try <- NULL
      }
      
      # If all PDE calls succeeded, compute morphology differences
      if (!inherits(PDEKernelsTry, "try-error") &&
          !inherits(pdeX1Try, "try-error") &&
          !inherits(pdeX2Try, "try-error")) {
        pdeX1 <- pdeX1Try$paretoDensity
        pdeX2 <- pdeX2Try$paretoDensity
        pdeDiff <- abs(matrixStats::rowDiffs(cbind(pdeX2, pdeX1)))
        
        n1 <- sum(Cls == UniqueCls[1])
        n2 <- sum(Cls == UniqueCls[2])
        Momentum1 <- sum(sign(Data[Cls == UniqueCls[1]]) *
                           log10(abs(Data[Cls == UniqueCls[1]]) + 1)) / log10(n1 + 1)
        Momentum2 <- sum(sign(Data[Cls == UniqueCls[2]]) *
                           log10(abs(Data[Cls == UniqueCls[2]]) + 1)) / log10(n2 + 1)
        DirMorph <- ifelse(Momentum2 < Momentum1, -1, 1)
        
        if (length(PDEKernels) == length(pdeDiff)) {
          MorphDiff <- caTools::trapz(PDEKernels, pdeDiff)
        }
        
        # Plot with PDE densities if requested
        if (PlotIt) {
          plot2Densities(Data, Cls,
                         PDEKernels = PDEKernels,
                         pdeX1 = pdeX1,
                         pdeX2 = pdeX2,
                         pde = pde,
                         meanLines = meanLines,
                         medianLines = medianLines,
                         col = col,
                         ...)
        }
      } else {
        # PDE failed, fallback to KDE plot if requested
        if (PlotIt) {
          plot2Densities(Data, Cls,
                         pde = FALSE,
                         meanLines = meanLines,
                         medianLines = medianLines,
                         col = col,
                         ...)
        }
      }
    }
    
    # Combine the weighted morphological and CT differences
    ImpactX2X1 <- CTDiffWeight * (DirCT * abs(CTDiff)) +
      (1 - CTDiffWeight) * (DirMorph * abs(MorphDiff))
  }
  
  # If KS-test suggests no significant difference, set impact to zero
  if (KSpval >= 0.05) ImpactX2X1 <- 0
  
  # Prepare density data frame for return
  if (is.null(PDEKernels) || is.null(pdeX1) || is.null(pdeX2)) {
    density_df <- data.frame(PDEKernels = numeric(0),
                             pde_Cls1 = numeric(0),
                             pde_Cls2 = numeric(0))
  } else {
    density_df <- cbind.data.frame(PDEKernels = PDEKernels,
                                   pde_Cls1 = pdeX1,
                                   pde_Cls2 = pdeX2)
  }
  
  list("Impact" = ImpactX2X1,
       "MorphDiff" = MorphDiff,
       "CTDiff" = CTDiff,
       "density_df" = density_df)
}
