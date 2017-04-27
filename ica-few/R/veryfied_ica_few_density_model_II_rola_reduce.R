rm(list = ls())
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE_II_grad.R")


m=c(0,0,0)
s=c(1,1,1)
t=c(1,1,1)
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
ica_kk=2
t=t[1:ica_kk]

#one<-diag(c(1,1,1))
param<-c(m,covariance)
#param<-c(m,e$vectors,s,t)

GSGaussian_par_mle_new(data=X, param=param, ica_k=ica_kk)

round(grad(func=ica_few_II, x=param, data=X, ica_k=ica_kk),digits=7)
round(gr_ica_few_II(param=param, data=X, ica_k=ica_kk),digits=7)

