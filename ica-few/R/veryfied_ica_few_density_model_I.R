rm(list = ls())
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")

fr_inv <- function(param, x) {  
  n=length(x[,1]); 
  d=length(x[1,]);  
  m<-matrix(param[1:d]); 
  w=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE);
  w=w+diag(x = 0.000000000001, nrow=d, ncol=d)
  vv=w;
  ans=1; 
  gAll<-c();
  s1All<-c();
  s2All<-c(); 
  tt<-c();
  ss<-c();  
  for(j in (1:d) ){ 
    xMove<-t(apply(x, 1, function(x) (vv[j,]%*%(matrix(x)-m))[1,1] )); 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3)); 
    ans=ans*g;
      
    s1All<-c(s1All,s1); 
    s2All<-c(s2All,s2);
    gAll<-c(gAll,g);
    lTemp=( (s1)^(2/3)*g/n); 
    lTemp=sqrt(lTemp); 
    tTemp=((s2/s1)^(1/3)); 
    tt<-c(tt,tTemp*tTemp); 
    ss<-c(ss,lTemp);
  }; 
  if(det(w)<=0){
    w<-w+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  
  #ans=ans*abs(abs(det(w))^(-2/3)); 
  
  ans=(d*n/2)*log((2*n)/(pi*exp(1))) -(3*n/2)*log(ans*abs(det(vv)^(-2/3)))  
  
  if(ans<=0){
    return(log(0.000000000001));
  }else{
    return(ans);
  }  
} 


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
  return(-ret- n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
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
#MLE_GSGaussian_par(data=X, param=param, ica_k=ica_kk)
#MLE_GSGaussian_par_mle(data=X, param=param)
#MLE_GSGaussian_par_mle_new(data=X, param=param, ica_k=ica_kk)

print(GSGaussian_par(data=X, param=param, ica_k=ica_kk))

#param<-c(m,covariance)
#fr_inv(param,x=X)
