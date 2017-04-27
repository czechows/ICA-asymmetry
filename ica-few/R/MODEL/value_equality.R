
rm(list = ls())
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
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=d,byrow=F);
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
#  ans2 = - ans2 - n*(1/2*log(det(cov(xMove))) )
  ans2 = - ans2 - n*(log(1/abs(det(vv))))
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



m=c(0,0)
s=c(1,1)
t=c(0.5,0.2)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=1000)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=1000)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)), 2, 2)
X <- S %*% A


m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_k=2

param<-c(m,e$vectors[,1:ica_k],s,t)
##########################################
MLE_GSGaussian_par(data=X, param=param, ica_k)
MLE_GSGaussian_par_mle(data=X, param=param)
##########################################


#########################################

m=c(0,0,0)
s=c(1,1,1)
t=c(0.5,0.2,1)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=1000)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=1000)
vawe3_l = rnorm(n=1000,mean=m[3],sd=1)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(4), 1/sqrt(4), -1/sqrt(4), 1/sqrt(4)), 2, 2)
#A <- matrix(c(1/sqrt(3), 1/sqrt(3), -1/sqrt(3), 1/sqrt(3)), 2, 2)
X <- S %*% A
X <- cbind(X,vawe3_l)
library("rgl")
plot3d(X,aps=1)


m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_k=3
s=s[1:ica_k]
t=t[1:ica_k]

param<-c(m,e$vectors[,1:ica_k],s,t)

##########################################
MLE_GSGaussian_par(data=X, param=param, ica_k)
MLE_GSGaussian_par_mle(data=X, param=param)
##########################################

s=c(1,1,1)
t=c(0.5,0.2,1)
ica_k=2
s=s[1:ica_k]
t=t[1:ica_k]
param<-c(m,e$vectors[,1:ica_k],s,t)
#MLE_GSGaussian_par(data=X, param=param, ica_k)
