\name{Impact}
\alias{Impact}
\title{Impact effect size measure}
\description{
  Calculates the Impact effect size measure that is based on the group distance and the difference in the shape of the data distribution between two groups.
}
\usage{
Impact(Data, Cls, PlotIt = FALSE, pde = TRUE, 
  col = c("red","blue"), meanLines = FALSE, medianLines = FALSE, ...)
}
\arguments{
 \item{Data}{Numeric vector containing the values of both groups.}
 \item{Cls}{Grouping vector or factor of same length as Data, containing exactly two distinct classes.}
 \item{PlotIt}{Logical; if TRUE, plots the probability density function (PDF) of the two groups using Pareto density estimation or kernel density as fallback.}
 \item{pde}{Logical; if TRUE, attempts Pareto Density Estimation (PDE) for the PDFs in the plot. Ignored if \code{PlotIt = FALSE}.}
 \item{col}{Character vector of length two specifying colors for the two groups in the plot. Ignored if \code{PlotIt = FALSE}.}
 \item{meanLines}{Logical; if TRUE, draws vertical lines at group means in the plot. Ignored if \code{PlotIt = FALSE}.}
 \item{medianLines}{Logical; if TRUE, draws vertical lines at group medians in the plot. Ignored if \code{PlotIt = FALSE}.}
 \item{\dots}{Further graphical parameters passed to the plotting function if \code{PlotIt = TRUE}.}
}
\value{
  Returns a list with the following components:
  \describe{
    \item{Impact}{Numeric scalar; the combined effect size measure based on difference in medians and distribution shapes.}
    \item{MorphDiff}{Numeric scalar; the extent of difference in shapes of the probability density functions.}
    \item{CTDiff}{Numeric scalar; the extent of difference in group medians.}
    \item{density_df}{Data frame with columns \code{PDEKernels}, \code{pde_Cls1}, \code{pde_Cls2} containing the density kernel points and Pareto densities for each group. This may be empty if Pareto density estimation was not successful or disabled.}
  }
}
\details{
  The Impact effect size measure combines central tendency differences (based on group medians) and morphological differences (based on Pareto density estimation).  
  If Pareto density estimation fails or is disabled, the density plots fallback to standard kernel density estimates.  
  The function can optionally plot these densities along with mean and/or median reference lines.
}
\references{
  Lotsch, J., and Ultsch, A. (2019): ImpactEffectsize -- an R Package for calculation and visualisation of the Impact distance and shape based effect size measure.
}
\author{
  Jorn Lotsch and Alfred Ultsch
}
\examples{
## Example 1: Use Impact with plotting
data("FeatureselectionData")
ImpactSize <- Impact(Data = FeatureselectionData$Var0011, 
  Cls = FeatureselectionData$Classes, PlotIt = TRUE)
  
## Example 2: Impact without plotting
ImpactSize <- Impact(Data = FeatureselectionData$Var0011, 
  Cls = FeatureselectionData$Classes, PlotIt = FALSE)

## example 2
data("BcellLymphomaCD79")
data("FeatureselectionData")
data("FlowcytometricData")
data("SameMeansData")
data("StocksFluctuation")
}
\keyword{Impact}
\keyword{effect size}
