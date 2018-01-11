#' LinearBaseline
#'
#'
#'
#' @name LinearBaseline
#' @author Jan Walach
#' @param x data matrix with rows as samples and columns as variables
#' @param prd.td cutoff parameter

#' @description Linear Baseline normalization.
#' @author Jan Walach
#' @details Linear Baseline normalization as described in [1].
#' @return The function returns object of type "matrix".
#' \item{x}{  matrix after normalization}

#' @export
#' @references [1] CAC Vol 82: Data Analysis for Omic Sciences: Methods and Applications chapters are due, Chapter Data normalization and scaling:  consequences forthe analysis in omics sciences
#' @examples
#' data(mcad)
#' x <- mcad[,-1]
#' LinearBaseline(x)

LinearBaseline <- function(x)
{
  x <- t(x)
  linear.baseline <- apply(x,1,median)
  baseline.mean <- mean(linear.baseline)
  sample.means <- apply(x,2,mean)
  linear.scaling <- baseline.mean/sample.means
  return( t(t(t(x)*linear.scaling)) )
}
