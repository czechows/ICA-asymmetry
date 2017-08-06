#' The ICA function based on asymmetry
#'
#' @return S, M, W, c as in the ICAA paper of Spurek et al. (2017)
#' @keywords Independent Component Analysis
#' @export
#' @examples 
#' ICAA(data)
ICAA<-function(X, noise=0, generalized=0, M = colMeans(X), W = eigen(cov(X),symmetric=TRUE)[[2]], c=2.0, minimum = -5.0, accuracy=1e-7 )
{
  # returns minimum, writes argmin to M, W, c
  ICA( X, M, W, c, noise, generalized, minimum, accuracy )

  return( list( scale(X%*%W), M, W, c) )
}


