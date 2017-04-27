rm(list = ls())
library("numDeriv")
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")

GSGaussian_par <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F); 
  ans=1
  xMove1<-t(apply(data, 1, function(x) ( (solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m)) ) );
  for(j in (1:ica_k) ){ 
    xMove=xMove1[,j] 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    ans=ans*g;
  } 
  ret=(ica_k*n/2)*log((2*n)/(pi*exp(1))) -(3*n/2)*log(ans)    
  #return(-ret- n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
  return(ret)
}

gr_GSGaussian_par <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  ans1=rep(0,d);   
  ans_c<-c()
  vv=(solve(t(v)%*%v)%*%t(v))
  xMove1<-t(apply(data, 1, function(x) ( vv%*%(matrix(x)-m)) ) );  
  for(j in (1:ica_k) ){ 
    xMove=xMove1[,j] 
    #xMove_N<-t(apply(data, 1, function(x) ( ((solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m))[j,]*v[j,] ) ) );     
    #xMove_N_w<-t(apply(data, 1, function(x) ( (2*(solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m))[j,]*(matrix(x)-m) ) ) );    
    xMove_N<-t(apply(data, 1, function(x) ((v[j,]%*%(matrix(x)-m))[1,1]*v[j,]) ));     
    xMove_N_w<-t(apply(data, 1, function(x) (2*(v[j,]%*%(matrix(x)-m))[1,1]*(matrix(x)-m) ) )); 
    
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    s1_N<-apply(xMove_N[xMove<=0,], 2, sum) 
    s2_N<-apply(xMove_N[xMove>0,], 2, sum) 
    s1_N_w<-apply(xMove_N_w[xMove<=0,], 2, sum) 
    s2_N_w<-apply(xMove_N_w[xMove>0,], 2, sum)
    
    g=(s1^(1/3)+s2^(1/3));     
    gg=(2/3)*(s1^(-2/3))*s1_N+(2/3)*(s2^(-2/3))*s2_N
    ggg=(1/3)*(s1^(-2/3))*s1_N_w+(1/3)*(s2^(-2/3))*s2_N_w 
    ans1=ans1-1/g*gg;
    ans_c<-c(ans_c,1/g*ggg)
  }
  ans=c(ans1, ans_c+ c((-2/3)))
  return(ans);
}

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
ica_kk=2
s=s[1:ica_kk]
t=t[1:ica_kk]

#param<-c(m,e$vectors[,1:ica_kk],s,t)
param<-c(m,covariance)
##########################################
print(GSGaussian_par(data=X, param=param, ica_k=ica_kk))

round(grad(func=GSGaussian_par, x=param, data=X, ica_k=ica_kk),digits=7)
round(gr_GSGaussian_par(param=param, data=X, ica_k=ica_kk),digits=7)


