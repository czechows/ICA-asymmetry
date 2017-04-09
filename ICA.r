Rcpp::sourceCpp('ICA.cpp')

myX = matrix( 1:6, nrow=3, ncol=2 )
myM = c(2,4)
myW = matrix( 1:4, nrow=2, ncol=2 )

ICA( myX, myM, myW )
