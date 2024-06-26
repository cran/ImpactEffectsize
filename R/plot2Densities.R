# Plots the pdf of two groups.
# Default is the Pareto density estimation (PDE).
# A standard pdf can be chosen or will be chosen in the case of zero variance of the data.
#' @export
#' @importFrom graphics abline lines plot
#' @importFrom stats density
#' @importFrom methods hasArg
#' @export
plot2Densities <- function(Data, Cls, col = c("red", "blue"), pde = TRUE, meanLines = FALSE, medianLines = FALSE, ...) {
  # Check if Data and Cls have the same length
  if (length(Data) != length(Cls)) {
    stop("Impact: Data and Cls have different lengths!")
  }
  
  # Get the unique classes
  UniqueCls <- sort(unique(Cls))
  
  # Check if there are enough data points for each class
  if (length(table(Data = Data[Cls == UniqueCls[1]])) < 2 ||
      length(table(Data = Data[Cls == UniqueCls[2]])) < 2) {
    suppressWarnings(pdeX1Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[1]]), TRUE))
    suppressWarnings(pdeX2Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[2]]), TRUE))
  } else {
    suppressWarnings(PDEKernelsTry <- try(ParetoDensityEstimationIE(Data), TRUE))
    if (!inherits(PDEKernelsTry, "try-error")) {
      PDEKernels <- PDEKernelsTry$kernels
      suppressWarnings(pdeX1Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[1]], kernels = PDEKernels), TRUE))
      suppressWarnings(pdeX2Try <- try(ParetoDensityEstimationIE(Data = Data[Cls == UniqueCls[2]], kernels = PDEKernels), TRUE))
    } else {
      message("Pareto density estimation failed. Reverting to standard pdf.")
      pde <- FALSE
    }
  }
  
  # Handle errors in Pareto density estimation
  if (inherits(pdeX1Try, "try-error")) {
    message("Pareto density estimation failed. Reverting to standard pdf.")
    pde <- FALSE
  }
  if (inherits(pdeX2Try, "try-error")) {
    message("Pareto density estimation failed. Reverting to standard pdf.")
    pde <- FALSE
  }
  
  # Calculate the density values
  if (pde == FALSE) {
    pdx1 <- stats::density(Data[Cls == UniqueCls[1]])$x
    pdx2 <- stats::density(Data[Cls == UniqueCls[2]])$x
    pd1 <- stats::density(Data[Cls == UniqueCls[1]])$y
    pd2 <- stats::density(Data[Cls == UniqueCls[2]])$y
  } else {
    pdx1 <- pdeX1Try$kernels
    pdx2 <- pdeX2Try$kernels
    pd1 <- pdeX1Try$paretoDensity
    pd2 <- pdeX2Try$paretoDensity
  }
  
  # Determine the plot limits
  xmin <- min(pdx1, pdx2)
  xmax <- max(pdx1, pdx2)
  ymax <- max(pd1, pd2)
  
  # Plot the densities
  graphics::plot(pd1 ~ pdx1, type = "l", lwd = 3, col = col[1], xlim = c(xmin, xmax), ylim = c(0, ymax), ...)
  graphics::lines(pd2 ~ pdx2, lwd = 3, col = col[2], ...)
  
  # Add median and mean lines if requested
  if (hasArg("medianLines") == TRUE && medianLines == TRUE) {
    graphics::abline(v = c_median(Data[Cls == UniqueCls[1]]), col = "magenta")
    graphics::abline(v = c_median(Data[Cls == UniqueCls[2]]), col = "magenta", lty = 2)
  }
  if (hasArg("meanLines") == TRUE && meanLines == TRUE) {
    graphics::abline(v = mean(Data[Cls == UniqueCls[1]]), col = "darkgreen")
    graphics::abline(v = mean(Data[Cls == UniqueCls[2]]), col = "darkgreen", lty = 2)
  }
}
