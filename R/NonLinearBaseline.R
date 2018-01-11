#' NonLinearBaseline
#'
#'
#'
#' @name NonLinearBaseline
#' @author Jan Walach
#' @param x data matrix with rows as samples and columns as variables
#' @param prd.td cutoff parameter

#' @description Non-Linear Baseline normalization.
#' @author Jan Walach
#' @details Non-Linear Baseline normalization as described in [1]. The function is based on function from Affy package [2]
#' @return The function returns object of type "matrix".
#' \item{liwong.x}{  matrix after normalization}

#' @export
#' @references [1] CAC Vol 82: Data Analysis for Omic Sciences: Methods and Applications chapters are due, Chapter Data normalization and scaling:  consequences forthe analysis in omics sciences
#' @references [2] L.Gautier and L. Cope and B.M. Bolstad and R. A. Irizarry (2004) affy---analysis of Affymetrix GeneChip data at the probe level, Bioinformatics, 20(3), 307--315, 2004
#' @examples
#' data(mcad)
#' x <- mcad[,-1]
#' NonLinearBaseline(x,prd.td=c(0.003,0.007))

NonLinearBaseline <- function(x,prd.td=c(0.003,0.007))
{
  x <- t(x)
  library(affy)

  # First step: Find baseline sample
  average.intensity <- apply(x,2,mean)
  median.number <- round(ncol(x)/2 + 0.1) # additional 0.1 ensures that it rounds properly
  ordering <- order(average.intensity)
  median.sample.number <- ordering[median.number]
  median.sample <- x[,median.sample.number]

  # Apply normalization
  liwong.x=vector()
  for(i in 1:ncol(x)){
    liwong.model <- normalize.invariantset(data=x[,i],
                                           ref=median.sample,
                                           prd.td=prd.td)
    liwong.sample <- predict(liwong.model$n.curve$fit,x[,i])
    liwong.x <- cbind(liwong.x,liwong.sample$y)
  }
  return(t(liwong.x))
}
