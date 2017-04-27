rm(list = ls())
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")

GSGaussian_par_mle_new <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  #v=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE)
  #p=d
  #q=ica_k
  ans=1
  ans1=1
  #vv=(solve(t(v)%*%v)%*%t(v))
  vv=v
  xMove1<-t(apply(data, 1, function(x) ( t(vv)%*%(matrix(x)-m)) ) );  
  for(j in (1:d) ){ 
    xMove=xMove1[,j] 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    if(j<=ica_k){
      ans=ans*g;
    }else{
      ans1=(ans1*(s1+s2)/n);
    }
    
  } 
  #param_temp<-c(m,v,ss,tt)
  #cc<-(n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1))
  #ret=cc+(-3*n/2)*(log(ans1^(1/3))+log(ans)+log(abs(det(v))^(2/3)))
  ret=( (ans1^(1/3))*(ans)*(abs(det(vv))^(-2/3))  )
  #ret=(log(ans1^(1/3))+log(ans)+log(abs(det(v))^(2/3)))
  if(ret<=0){
    return(log(0.000000000001));
  }else{
    return(log(ret))
  }
}



gr_GSGaussian_par_mle_new <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  #v=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE)
  #vv=(solve(t(v)%*%v)%*%t(v))
  vv=v
  xMove1<-t(apply(data, 1, function(x) (t(vv)%*%(matrix(x)-m))));
  ans1=rep(0,d);   
  ans_c<-c()
  for(j in (1:d) ){ 
    xMove=xMove1[,j] 
    xMove_N<-t(apply(data, 1, function(x) ( (t(vv)%*%(matrix(x)-m))[j,]*t(vv)[j,] ) ) );     
    xMove_N_w<-t(apply(data, 1, function(x) ( (2*t(vv)%*%(matrix(x)-m))[j,]*(matrix(x)-m) ) ) );  
    #xMove_N_w<-t(apply(data, 1, function(x) (2*(vv[j,]%*%(matrix(x)-m))[1,1]*(matrix(x)-m) ) ));  
    
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    s1_N<-apply(xMove_N[xMove<=0,], 2, sum) 
    s2_N<-apply(xMove_N[xMove>0,], 2, sum) 
    s1_N_w<-apply(xMove_N_w[xMove<=0,], 2, sum) 
    s2_N_w<-apply(xMove_N_w[xMove>0,], 2, sum)
    
    g=(s1^(1/3)+s2^(1/3));     
    gg=(2/3)*(s1^(-2/3))*s1_N+(2/3)*(s2^(-2/3))*s2_N
    gg_1=(2/(3*(s1+s2)))*(s1_N+s2_N)
    ggg=(1/3)*(s1^(-2/3))*s1_N_w+(1/3)*(s2^(-2/3))*s2_N_w
    ggg_1=(1/(3*(s1+s2)))*(s1_N_w+s2_N_w) 
    if(j<=ica_k){
      ans1=ans1-1/g*gg;
      ans_c<-c(ans_c,1/g*ggg)
    }else{
      ans1=ans1-gg_1;
      ans_c<-c(ans_c,ggg_1)
    }    
    
  }
  ans=c(ans1, ans_c+ c((-2/3)*(solve(t(vv)))) )
  return(ans);
}

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
ica_kk=3
t=t[1:ica_kk]

#one<-diag(c(1,1,1))
param<-c(m,covariance)
#param<-c(m,e$vectors,s,t)

GSGaussian_par_mle_new(data=X, param=param, ica_k=ica_kk)

round(grad(func=GSGaussian_par_mle_new, x=param, data=X, ica_k=ica_kk),digits=7)
round(gr_GSGaussian_par_mle_new(param=param, data=X, ica_k=ica_kk),digits=7)




