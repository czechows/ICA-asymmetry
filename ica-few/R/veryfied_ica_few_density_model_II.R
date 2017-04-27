#rm(list = ls())
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")

GSGaussian_par_mle_new <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  
  s=matrix(param[(d+d^2+1):(2*d+d^2)]);
  t=matrix(param[(2*d+d^2+1):(2*d+ica_k+d^2)]);
  p=d
  q=ica_k
  t=c(t,rep(1,d-ica_k))
  ans=1
  ans1=1
  xMove1<-t(apply(data, 1, function(x) ( (solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m)) ) );
  print("...")
  gAll<-c();
  s1All<-c();
  s2All<-c(); 
  tt<-c();
  ss<-c();   
  for(j in (1:d) ){ 
    xMove=xMove1[,j] 
    xMove<-t(apply(data, 1, function(x) (v[,j]%*%(matrix(x)-m))[1,1] ));
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    
    if(j<=ica_k){
      tTemp=((s2/s1)^(1/3));    
      lTemp=( (s1)^(2/3)*g/n);
      lTemp=sqrt(lTemp);
      tt<-c(tt,tTemp);
      ss<-c(ss,lTemp);
      ans=ans*g;
    }else{
      lTemp=((s1+s2)/n)
      lTemp=sqrt(lTemp); 
      ss<-c(ss,lTemp);
      ans1=(ans1*(s1+s2)/n);
    }
  } 
  param_temp<-c(m,v,ss,tt)
  print("...")
  print(MLE_GSGaussian_par_mle_new(data=data, param=param_temp, ica_k=ica_k))
  t=c(tt,1)
  param<-c(m,v,ss,t)
  print(MLE_GSGaussian_par_mle(data=data, param=param))
  print("...")
  print("_____")
  MLE<-1
  print(MLE)
  print("///")
  mle_pre=log( (ans1^(1/3))*(ans)*(abs(det(v))^(-2/3)) )
  print( mle_pre )
  print("_____")
  ret=(n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1)) -(3*n/2)*log(ans)-(3*n/2)*log(ans1^(1/3)) 
  print( (n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1)) 
         -(3*n/2)*  (  log(ans*ans1^(1/3)*(abs(det(v))^(-2/3)))  ) 
         )
  print( (n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1)) 
         -(3*n/2)*  (  mle_pre  ) 
  )
  #ret=ret+log(abs(abs(det(v))))
  #print(-ret - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
  return(-ret - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
}

m=c(0,0,0)
s=c(1,1,1)
t=c(1,1,1)
vawe1_l = rspliNorm(m=m[1],s=s[1],t=t[1],size=2000)
vawe2_l = rspliNorm(m=m[2],s=s[2],t=t[2],size=2000)
vawe3_l = rnorm(n=2000,mean=m[3],sd=1)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1, 2, -1, 1), 2, 2)
X <- S %*% A
X <- cbind(X,vawe3_l)
plot3d(X,aps=1)


m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_kk=3
t=t[1:ica_kk]

#covariance<-matrix(c(4.9748, 1.0052, 0.0754, 1.0052, 2.0103, 0.0530, 0.0754, 0.0530, 1.0284),3)
#one<-diag(c(1,1,1))
param<-c(m,e$vectors,s,t)
#param<-c(m,e$vectors,s,t)

#ans2=ans2+log(abs(abs(det(w))^(-2/3)))
GSGaussian_par_mle_new(data=X, param=param, ica_k=ica_kk)
pppp(data=X, param=param, ica_k=ica_kk)

library(Rmixmod)
mixmodCluster(data.frame(X),1)
 