Rcpp::sourceCpp('ICA.cpp')

myX = matrix( runif(n=16, min=-1, max=1), nrow=8, ncol=2 )

myM = runif( n=2, min=-1, max=1 ) 
myW = matrix( runif(n=4,min=-1,max=1), nrow=2, ncol=2 )

ICA( myX, myM, myW )
