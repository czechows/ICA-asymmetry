#' The ICA function based on asymmetry
#'
#' @return S, M and W as in the ICAA paper of Spurek et al. (2017)
#' @keywords Independent Component Analysis
#' @export
#' @examples 
#' ICAA(data)
ICAA<-function(X)
{
  M = colMeans(X) 
  W = eigen(cov(X),symmetric=TRUE)[[2]]
  
  # returns minimum, writes argmin to M, W
  ICA( X, M, W )
  
  return( list( scale(X%*%W), M, W) )
}


