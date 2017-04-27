
rm(list = ls())
library("rgl")
library("numDeriv")
set.seed(1234212);
norm_vec <- function(x) sqrt(sum(x^2))

################################################
# funkcja gestosci GSG
################################################
SGaussian1D<- function(x, m, l, t){
  c<-sqrt(2/pi)*(1/l)*(1/(1+t))
  if(x<=m){
    ans=c*exp(-(1/(2*l^2))*(x-m)^2)    
  }else{
    ans=c*exp(-(1/(2*t^2*l^2))*(x-m)^2)
  }
  return(ans)
}

GSGaussian<- function(x, v, m, s, t){
  ans<-1;
  d=length(x); 
  if(det(v)==0){
    v<-v+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  vv=solve(v)
  #
  for(j in (1:length(x)) ){
    ans=ans*SGaussian1D(x=vv[j,]%*%(matrix(as.numeric(x))-m),m=0,l=(s[j]),t=t[j])
    #    print(SGaussian1D(x=vv[j,]%*%(matrix(as.numeric(x))-m),m=0,l=(s[j]),t=t[j]))
  }
  return(1/abs(det(v))*ans)
}

GSGaussian_ort<- function(x, v, m, s, t){
  ans<-1;
  #
  for(j in (1:length(x)) ){
    ans=ans*SGaussian1D(x=x[j],m=0,l=(s[j]),t=t[j])
  }
  return(ans)
}

MLE_GSGaussian <- function(x, v, m, s, t){
  ans2=0
  for(i in (1:length(x[,1])) ){
    el=x[i,]
    #    ans2=ans2-log(logistic(x=el, v=v, m=m, s=s))
    ans2=ans2+log(GSGaussian(x=el, v=v, m=m, s=s, t=t))
  }
  return(ans2)
}

MLE_GSGaussian_par_mle <- function(data, param){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  s=matrix(param[(d+d^2+1):(2*d+d^2)]);
  t=matrix(param[(2*d+d^2+1):(3*d+d^2)]);
  ans2=0
  for(i in (1:length(data[,1])) ){
    el=data[i,]
    #    ans2=ans2-log(logistic(x=el, v=v, m=m, s=s))
    if(GSGaussian(x=el, m=m, v=v, s=s, t=t)<=0){
      ans2=ans2+log(0.0000000001)
    }else{
      ans2=ans2+log(GSGaussian(x=el, m=m, v=v, s=s, t=t))    
    }
    #ans2=ans2+log(GSGaussian(x=el, v=v, m=m, s=s, t=t))
  }
  
  return(-ans2)
}

MLE_GSGaussian_par_mle_new <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  s=matrix(param[(d+d^2+1):(2*d+d^2)]);
  t=matrix(param[(2*d+d^2+1):(2*d+ica_k+d^2)]);
  t=c(t,rep(1,d-ica_k))
  ans2=0
  for(i in (1:length(data[,1])) ){
    el=data[i,]
    #    ans2=ans2-log(logistic(x=el, v=v, m=m, s=s))
    if(GSGaussian(x=el, m=m, v=v, s=s, t=t)<=0){
      ans2=ans2+log(0.0000000001)
    }else{
      ans2=ans2+log(GSGaussian(x=el, m=m, v=v, s=s, t=t))    
    }
    #ans2=ans2+log(GSGaussian(x=el, v=v, m=m, s=s, t=t))
  }
  
  return(-ans2)
}


MLE_GSGaussian_par <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  s=matrix(param[(d+d*ica_k+1):(d+ica_k+d*ica_k)]);
  t=matrix(param[(d+ica_k+d*ica_k+1):(d+2*ica_k+d*ica_k)]);
  
  
  ans2=0
  vv=v#solve(v)
  xMove<-t(apply(data, 1, function(x) ( (solve(t(vv)%*%vv)%*%t(vv))%*%(matrix(x)-m)) ) );
  for(i in (1:length(xMove[,1])) ){
    el=xMove[i,]
    #    ans2=ans2-log(logistic(x=el, v=v, m=m, s=s)) 
    if(GSGaussian_ort(x=el, m=m, v=vv, s=s, t=t)<=0){
      ans2=ans2+log(0.0000000001)
    }else{
      ans2=ans2+log(GSGaussian_ort(x=el, m=m, v=vv, s=s, t=t))    
    }
  }  
  #ans2 = - ans2 - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove))) )
  ans2 = - ans2 - n*(1/2*log(det(cov(xMove))) )
  #ans2 = - ans2 - n*(log(1/abs(det(vv))))
  #ans2 = - ans2 - n*(log(det(solve((vv)%*%t(vv))%*%(vv))) ) 
  return(ans2)
}

###############################################################################################
###############################################################################################

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


one<-diag(c(1,1,1))



m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=200)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=200)
vawe3_l = rnorm(n=200,mean=m[3],sd=1)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(4), 1/sqrt(4), -1/sqrt(4), 1/sqrt(4)), 2, 2)
#A <- matrix(c(1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)), 2, 2)
X <- S %*% A
X <- cbind(X,vawe3_l)
plot3d(X,aps=1)


m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_k=3
s=s[1:ica_k]
t=t[1:ica_k]

param<-c(m,e$vectors[,1:ica_k],s,t)
#param<-c(m,one[,1:ica_k],s,t)

##########################################
#MLE_GSGaussian_par(data=X, param=param, ica_k)
#MLE_GSGaussian_par_mle(data=X, param=param)
#MLE_GSGaussian_par_mle_new(data=X, param=param, ica_k)
##########################################


op<-optim(par=param, fn=MLE_GSGaussian_par, data=X, ica_k=ica_k)#, method="SANN")
op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
s=matrix(op$par[(d+d*ica_k+1):(d+ica_k+d*ica_k)]);
t=matrix(op$par[(d+ica_k+d*ica_k+1):(d+2*ica_k+d*ica_k)]);

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


m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
t=t[1:ica_k]

param<-c(m,e$vectors,s,t)
#param<-c(m,one,s,t)
##########################################
#MLE_GSGaussian_par(data=X, param=param)
#MLE_GSGaussian_par_grad(x=X, param=param)
#MLE_GSGaussian_par_mle(x=X, param=param)
##########################################

op<-optim(param, MLE_GSGaussian_par_mle_new, data=X, ica_k=ica_k)
#op<-optim(par=param, fn=MLE_GSGaussian_par, gr=MLE_GSGaussian_par_grad, data=X)#, method="SANN")

op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d^2)],ncol=d,byrow=F);
s=matrix(op$par[(d+d^2+1):(2*d+d^2)]);
t=matrix(op$par[(2*d+d^2+1):(2*d+ica_k+d^2)]);
t=c(t,rep(1,d-ica_k))

v1<-v[1:3,1]
v2<-v[1:3,2]
v3<-v[1:3,3]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="blue", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="blue", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="blue", pch=20, lwd=2)

##########################################

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)

param<-c(m,e$vectors,s,t)
#param<-c(m,one,s,t)

op<-optim(param, MLE_GSGaussian_par_mle, data=X)

op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d^2)],ncol=d,byrow=F);
s=matrix(op$par[(d+d^2+1):(2*d+d^2)]);
t=matrix(op$par[(2*d+d^2+1):(3*d+d^2)]);

v1<-v[1:3,1]
v2<-v[1:3,2]
v3<-v[1:3,3]


segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="green", pch=20, lwd=6)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="green", pch=20, lwd=6)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="green", pch=20, lwd=6)
