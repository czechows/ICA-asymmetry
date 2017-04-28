#' The ICA function based on asymmetry
#'
#' @return M and W as in the ICAA paper
#' @keywords Independent Component Analysis
#' @export
#' @examples
#' ICA(data)
ICAA<-function(myX)
{
  myM = runif( n=3, min=-1, max=1 ) 
  myW = matrix( runif(n=9,min=-1,max=1), nrow=3, ncol=3 )
  
  # returns minimum, writes argmin to myM, myW
  ICA( myX, myM, myW )
  
  return(list(myM,myW))
}


