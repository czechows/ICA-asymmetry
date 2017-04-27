rm(list = ls())
set.seed(5)
#
#install.packages("matrixcalc")
#library(matrixcalc)
source("F:\\AfCEC\\artykuly\\ACEC\\R\\OSTATECZNY_MODEL\\GeneralSplitGausian.R")
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE_II_grad.R")

rspliNorm<-function(m,s,t,size){
  n <- 100000
  X <- runif(n,min=-10,max=10)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(n,min=0,max=1)
  ans<-c()
  for(i in (1:length(X)) ){
    if( Y[i] < SGaussian1D(X[i], m=m, l=s, t=t)  ){
      ans<-c(ans,X[i])
    }
  }
  return(ans[1:size])
}

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=2000)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=2000)
vawe3_l = rnorm(n=2000,mean=m[3],sd=1)


S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(4), 1/sqrt(4), -1/sqrt(4), 1/sqrt(4)), 2, 2)
X <- S %*% A
X <- cbind(X,vawe3_l)
plot3d(X,aps=1)


model_1<-convert_param_few_II(x=X, ica_k=2, n_starts=4) 

m<-model_1$m

v1<-model_1$cov[1,1:3]
v2<-model_1$cov[2,1:3]
v3<-model_1$cov[3,1:3]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="red", pch=20, lwd=2)