norm_vec <- function(x) sqrt(sum(x^2))

################################################
# G³owna funkcja
################################################
convert_param <- function(param, x) { # param - punkt startoey # x - dane  
  #op<-optim(param, fr_inv, gr=gr_fr, x=x,method="BFGS"); # minimalizacja grafjetowa
  op<-optim(param, fr_inv, x=x)
  n=length(x[,1])
  d=length(x[1,])  
  c1=sqrt(2/pi)^d
  m<-matrix(op$par[1:d])
  v=matrix(op$par[(d+1):length(op$par)],nrow=d,byrow=TRUE)
  if(det(v)<=0){
    v<-v+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  vv=solve(v)
  t<-c();
  s<-c();
  gAll<-c();
  s1All<-c();
  s2All<-c();        
  for(j in (1:d) ){ xMove<-t(apply(x, 1, function(x) (vv[j,]%*%(matrix(x)-m))[1,1] )); s1<-sum(xMove[xMove<=0]^2); s1=max(s1,0.000000000001); s2<-sum(xMove[xMove>0]^2); s2=max(s2,0.000000000001); s1All<-c(s1All,s1); s2All<-c(s2All,s2); g=(s1^(1/3)+s2^(1/3)); gAll<-c(gAll,g); lTemp=( (s1)^(2/3)*g/n); lTemp=sqrt(lTemp); tTemp=((s2/s1)^(1/3)); t<-c(t,tTemp); s<-c(s,lTemp);};
  if(op$value<=0 | is.na(op$value) | !is.finite(100)){
    op$value=0.000000000001
  }
  MLE=(n*d/2)*log((2*n)/(pi*exp(1))) + (-3*n/2)*( log( op$value ));    
  
  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=v,s=s,t=t,s1=s1All,s2=s2All,MLE=MLE,par=op$par,X=move_data))
}  
################################################
# funkcja, ktora minimalizujemy
################################################
fr_inv <- function(param, x) {  
  n=length(x[,1]); 
  d=length(x[1,]);  
  m<-matrix(param[1:d]); 
  w=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE);
  w=w+diag(x = 0.000000000001, nrow=d, ncol=d)
  vv=w;
  ans=1; 
  for(j in (1:ica_k) ){ 
    #xMove<-t(apply(x, 1, function(x) (vv[j,]%*%(matrix(x)-m))[1,1] )); 
    xMove<-t(apply(x, 1, function(x) ( solve((vv)%*%t(vv))%*%(vv%*%(matrix(x)-m)) ) ));
    xMove=xMove[,j]
    
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3)); 
    if(is.finite( ans*g )==T){
      ans=ans*g;
    }
    if(is.finite( ans )==F){
      print(j)
      #print(d)
      #print(s1)
      #print(s2)
      print(ans)
      print(log(ans*abs(abs(det(w))^(-2/3))))
      break;
    }
  }; 
  if(det(w)<=0){
    w<-w+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  
  ans=ans; 
  
  
  
  if(ans<=0){
    return(log(0.000000000001));
  }else{
    return(log(ans));
  }
  
} 
################################################
# pochodna tej funkcji
################################################
gr_fr <- function(param, x) {  
  n=length(x[,1]); 
  d=length(x[1,]);  
  m<-matrix(param[1:d]); 
  vv=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE);
  ans1=rep(0,d);   
  ans_c<-c()
  for(j in (1:d) ){ 
    xMove<-t(apply(x, 1, function(x) (vv[j,]%*%(matrix(x)-m))[1,1] ));
    xMove_N<-t(apply(x, 1, function(x) ((vv[j,]%*%(matrix(x)-m))[1,1] *vv[j,]) ));     
    xMove_N_w<-t(apply(x, 1, function(x) (2*(vv[j,]%*%(matrix(x)-m))[1,1]*(matrix(x)-m) ) ));    
    
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    s1_N<-apply(xMove_N[xMove<=0,], 2, sum)#sum(xMove_N[xMove<=0,]); 
    s2_N<-apply(xMove_N[xMove>0,], 2, sum)#sum(xMove_N[xMove>0,]); 
    s1_N_w<-apply(xMove_N_w[xMove<=0,], 2, sum)#sum(xMove_N_v[xMove<=0]); 
    s2_N_w<-apply(xMove_N_w[xMove>0,], 2, sum)#sum(xMove_N_v[xMove>0]);
    
    g=(s1^(1/3)+s2^(1/3));     
    gg=(2/3)*(s1^(-2/3))*s1_N+(2/3)*(s2^(-2/3))*s2_N
    ggg=(1/3)*(s1^(-2/3))*s1_N_w+(1/3)*(s2^(-2/3))*s2_N_w 
    ans1=ans1-1/g*gg;
    ans_c<-c(ans_c,1/g*ggg)
  }
  if(det(vv)<=0){
    vv<-vv+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  ans=c(ans1, ans_c+ c((-2/3)*(solve(vv))))
  return(ans);
}

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
  if(det(v)<=0){
    v<-v+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  vv=solve(v)
  #
  for(j in (1:length(x)) ){
    ans=ans*SGaussian1D(x=vv[j,]%*%(matrix(as.numeric(x))-m),m=0,l=(s[j]),t=t[j])
  }
  return(1/abs(det(v))*ans)
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
vawe1_l = rspliNorm(m=0,s=s[1],t=t[1],size=1000)
vawe2_l = rspliNorm(m=0,s=s[2],t=t[2],size=1000)

S <- cbind(scale(vawe1_l), scale(vawe2_l))
A <- matrix(c(1/sqrt(2), 1/sqrt(2), -1/sqrt(2), 1/sqrt(2)), 2, 2)
X <- S %*% A


vawe_sum = X[,1]
vawe_div = X[,2]

data1<-X[,1]
data2<-X[,2]

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
ica_k=2


param<-c(m,e$vectors[1:ica_k,])

model<-convert_param(param=param, x=X)  

covT<-model$cov
covT[,1]<-1/norm_vec(covT[,1])*covT[,1]
covT[,2]<-1/norm_vec(covT[,2])*covT[,2]
dataSGD<-scale(model$X%*%t(solve(covT))) 

#model_cov=(solve(covT))
#model_cov[,1]<-1/norm_vec(model_cov[,1])*model_cov[,1]
#model_cov[,2]<-1/norm_vec(model_cov[,2])*model_cov[,2]

model_cov=covT
model_cov[,1]<-1/norm_vec(model_cov[,1])*model_cov[,1]
model_cov[,2]<-1/norm_vec(model_cov[,2])*model_cov[,2]

#  model_cov_1=t(solve(covT))
#  model_cov_1[,1]<-1/norm_vec(model_cov_1[,1])*model_cov_1[,1]
#  model_cov_1[,2]<-1/norm_vec(model_cov_1[,2])*model_cov_1[,2]

plot(X, main = "",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)

