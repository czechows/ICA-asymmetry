rm(list = ls())
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=20)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=20)
vawe3_l = rnorm(n=20,mean=m[3],sd=1)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(4), 1/sqrt(4), -1/sqrt(4), 1/sqrt(4)), 2, 2)
X <- S %*% A
X <- cbind(X,vawe3_l)
plot3d(X,aps=1)


m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_kk=3
s=s[1:ica_kk]
t=t[1:ica_kk]

#one<-diag(c(1,1,1))
param<-c(m,e$vectors[,1:ica_kk],s,t)
#param<-c(m,one[,1:ica_k],s,t)

##########################################
#MLE_GSGaussian_par(data=X, param=param, ica_k=ica_kk)
#MLE_GSGaussian_par_mle(data=X, param=param)
#MLE_GSGaussian_par_mle_new(data=X, param=param, ica_k=ica_kk)

#MLE_GSGaussian_par_grad(data=X, param=param, ica_k=ica_kk)
#MLE_GSGaussian_par_mle_grad(data=X, param=param)
#MLE_GSGaussian_par_mle_new_grad(data=X, param=param, ica_k=ica_kk)

##########################################

op<-optim(par=param, fn=MLE_GSGaussian_par, data=X, ica_k=ica_kk)#, method="SANN")
#op<-optim(par=param, fn=MLE_GSGaussian_par, gr=MLE_GSGaussian_par_grad, data=X, ica_k=ica_kk, method="BFGS")#, method="SANN")
#op<-optim(par=param, fn=MLE_GSGaussian_par, gr=MLE_GSGaussian_par_mle_new_grad, data=X, ica_k=ica_kk, method="BFGS")#, method="SANN")


op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d*ica_kk)],ncol=ica_kk,byrow=F);
s=matrix(op$par[(d+d*ica_kk+1):(d+ica_kk+d*ica_kk)]);
t=matrix(op$par[(d+ica_kk+d*ica_kk+1):(d+2*ica_kk+d*ica_kk)]);

v1<-v[1:3,1]
v2<-v[1:3,2]
#v3<-v[1:3,3]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="red", pch=20, lwd=2)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="red", pch=20, lwd=2)

#segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
#           y=as.vector(t(c(m[2],m[2]+v3[2]))),
#           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="red", pch=20, lwd=4)
#segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
#segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)

###########################################################################################################

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
t=t[1:ica_kk]

param<-c(m,e$vectors,s,t)
#param<-c(m,one,s,t)
##########################################
#MLE_GSGaussian_par(data=X, param=param)
#MLE_GSGaussian_par_grad(x=X, param=param)
#MLE_GSGaussian_par_mle(x=X, param=param)
##########################################

op<-optim(param, MLE_GSGaussian_par_mle_new, data=X, ica_k=ica_kk)
#op<-optim(par=param, fn=MLE_GSGaussian_par_mle_new, gr=MLE_GSGaussian_par_mle_new_grad, data=X, ica_k=ica_kk, method="BFGS")

op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d^2)],ncol=d,byrow=F);
s=matrix(op$par[(d+d^2+1):(2*d+d^2)]);
t=matrix(op$par[(2*d+d^2+1):(2*d+ica_kk+d^2)]);
t=c(t,rep(1,d-ica_kk))

v1<-v[1:3,1]
v2<-v[1:3,2]
v3<-v[1:3,3]

norm_vec(v1)

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="green", pch=20, lwd=2)

##########################################
###########################################################################################################

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)


param<-c(m,e$vectors,s,t)
##########################################
#MLE_GSGaussian_par(data=X, param=param)
#MLE_GSGaussian_par_grad(x=X, param=param)
#MLE_GSGaussian_par_mle(x=X, param=param)
##########################################

op<-optim(param, MLE_GSGaussian_par_mle, data=X)
#op<-optim(par=param, fn=MLE_GSGaussian_par_mle, gr=MLE_GSGaussian_par_mle_grad, data=X, method="BFGS")


op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d^2)],ncol=d,byrow=F);
s=matrix(op$par[(d+d^2+1):(2*d+d^2)]);
t=matrix(op$par[(2*d+d^2+1):(3*d+d^2)]);
#v1<-v[1,1:3]
#v2<-v[2,1:3]
#v3<-v[3,1:3]

v1<-v[1:3,1]
v2<-v[1:3,2]
v3<-v[1:3,3]

norm_vec(v1)

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="blue", pch=20, lwd=6)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="blue", pch=20, lwd=6)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="blue", pch=20, lwd=6)






