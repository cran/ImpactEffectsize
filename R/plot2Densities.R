#Plots the pdf of two groups. 
#Default is the Pareto density estimation (PDE). 
#A standard pdf can be chosen or will be chosen in the case of zero variance of the data. 
#' @importFrom DataVisualizations ParetoDensityEstimation
#' @importFrom graphics abline lines plot
#' @importFrom stats density median
#' @importFrom methods hasArg
#' @export
plot2Densities <- function(Data,Cls,col=c("red","blue"), pde=TRUE, 
                           meanLines=FALSE,medianLines=FALSE,...){
  if(length(Data) != length(Cls)) stop("Impact: Data and Cls have different lengths!")
  if(var(Data[Cls==sort(unique(Cls))[1]]) == 0 | var(Data[Cls==sort(unique(Cls))[2]]) == 0) {
    message("One or both group data have no variance. Reverting to standard pdf")
    pde = FALSE
  }
  if(hasArg(pde) == TRUE & pde == FALSE) {
    pdx1 = density(Data[Cls==sort(unique(Cls))[1]])$x
    pdx2 = density(Data[Cls==sort(unique(Cls))[2]])$x
    pd1 = density(Data[Cls==sort(unique(Cls))[1]])$y
    pd2 = density(Data[Cls==sort(unique(Cls))[2]])$y
  } else {
    pdx1 <-DataVisualizations::ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[1]])$kernels
    pdx2 <- DataVisualizations::ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[2]])$kernels
    pd1 <-DataVisualizations::ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[1]])$paretoDensity
    pd2 <- DataVisualizations::ParetoDensityEstimation(Data = Data[Cls==sort(unique(Cls))[2]])$paretoDensity
  }
  xmin=min(pdx1,pdx2)
  xmax=max(pdx1,pdx2)
  ymax=max(pd1,pd2)
  plot(pd1~pdx1,type="l", lwd=3, col=col[1], xlim=c(xmin,xmax),ylim=c(0,ymax))
  lines(pd2~pdx2,lwd=3,col=col[2],...)
  if(hasArg(medianLines) == TRUE & medianLines == TRUE) {
    abline(v=median(Data[Cls==sort(unique(Cls))[1]]), col = "magenta")
    abline(v=median(Data[Cls==sort(unique(Cls))[2]]), col = "magenta",lty=2)
  }
  if(hasArg(meanLines) == TRUE & meanLines == TRUE) {
    abline(v=mean(Data[Cls==sort(unique(Cls))[1]]), col = "darkgreen")
    abline(v=mean(Data[Cls==sort(unique(Cls))[2]]), col = "darkgreen",lty=2)
  }
}
