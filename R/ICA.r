Rcpp::sourceCpp('src/ICA.cpp')

myX = matrix( runif(n=24, min=-1, max=1), nrow=8, ncol=3 )

myM = runif( n=3, min=-1, max=1 ) 
myW = matrix( runif(n=9,min=-1,max=1), nrow=3, ncol=3 )

# returns minimum, writes argmin to myM, myW
ICA( myX, myM, myW )

print( myM )
print( myW )



