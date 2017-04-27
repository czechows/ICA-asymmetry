
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
  vv=v#solve(v)
  #
  for(j in (1:length(x)) ){
    ans=ans*SGaussian1D(x=vv[j,]%*%(matrix(as.numeric(x))-m),m=0,l=(s[j]),t=t[j])
    #    print(SGaussian1D(x=vv[j,]%*%(matrix(as.numeric(x))-m),m=0,l=(s[j]),t=t[j]))
  }
  return(abs(det(v))*ans)
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


MLE_param_GeneralGaussian_1 <- function(x, v, m, s, t){
  ans2=0
  for(i in (1:length(x[,1])) ){
    el=x[i,]
    #    ans2=ans2-log(logistic(x=el, v=v, m=m, s=s))  
    ans2=ans2+log(GSGaussian(x=el, m=m, v=v, s=s, t=t))    
  }
  return(ans2)
}
###################
# funkcja tylko do sprawdzania poprawnosci wzorow
###################

GSGaussian_par_mle_new_veryfi <- function(data, param, ica_k){
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
  gAll<-c();
  s1All<-c();
  s2All<-c(); 
  tt<-c();
  ss<-c();   
  for(j in (1:d) ){ 
    xMove=xMove1[,j] 
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
      ans1=(ans1*sqrt((s1+s2)/n));
    }
    
  } 
  param_temp<-c(m,v,ss,tt)
  print("...")
  print(MLE_GSGaussian_par_mle_new(data=data, param=param_temp, ica_k=ica_k))
  t=c(tt,1)
  param<-c(m,e$vectors,ss,t)
  print(MLE_GSGaussian_par_mle(data=data, param=param))
  print("...")
  ret=(n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1)) -(3*n/2)*log(ans)-n*log(ans1) 
  print(-ret - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
  return(-ret - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
}

###################
# funkcja tylko do sprawdzania poprawnosci wzorow
###################
GSGaussian_par_veryfi <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  
  ans=1
  xMove1<-t(apply(data, 1, function(x) ( (solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m)) ) );
  
  gAll<-c();
  s1All<-c();
  s2All<-c(); 
  tt<-c();
  ss<-c();   
  for(j in (1:ica_k) ){ 
    xMove=xMove1[,j] 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    
    lTemp=( (s1)^(2/3)*g/n); 
    lTemp=sqrt(lTemp); 
    tTemp=((s2/s1)^(1/3));
    
    ss<-c(ss,lTemp);
    tt<-c(tt,tTemp);
    ans=ans*g;
  } 
  #m<-matrix(param[1:d]); 
  #v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  param_temp<-c(m,v,ss,tt)
  print("...")
  #print(MLE_GSGaussian_par_mle(data=X, param=param_temp))
  #print(-MLE_param_GeneralGaussian_1(x=X, v=v, m=m, s=ss, t=tt)) 
  print(MLE_GSGaussian_par(data=X, param=param_temp, ica_k=ica_k))
  #  MLE_GSGaussian_par_mle_new(data=X, param=param_temp, ica_k=ica_kk)
  print("...")
  ret=(ica_k*n/2)*log((2*n)/(pi*exp(1))) -(3*n/2)*log(ans)    
  return(-ret- n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
}



MLE_GSGaussian_par_mle_grad <- function(data, param){
  grad(func=MLE_GSGaussian_par_mle, x=param, data=data)
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
  #return(-ans2- n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(data))) ))
  return(-ans2)
}

MLE_GSGaussian_par_mle_new_grad <- function(data, param, ica_k){
  grad(func=MLE_GSGaussian_par_mle_new, x=param, data=data, ica_k=ica_k)
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
  #ans2=ans2+log(abs(det(v)))
  return(-ans2)
  #return(-ans2- n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(data)))))
}

#####################################################################
# funkcja liczaca ICA wersja I z wyliczonym wzorem
#####################################################################
GSGaussian_par_grad <- function(data, param, ica_k){
  grad(func=GSGaussian_par, x=param, data=data, ica_k=ica_k)
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

#####################################################################
# funkcja liczaca ICA wersja II z wyliczonym wzorem
#####################################################################
GSGaussian_par_mle_new_grad <- function(data, param, ica_k){
  grad(func=GSGaussian_par_mle_new, x=param, data=data, ica_k=ica_k)
}

GSGaussian_par_mle_new <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  p=d
  q=ica_k
  ans=1
  ans1=1
  xMove1<-t(apply(data, 1, function(x) ( (solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m)) ) );
  gAll<-c();
  s1All<-c();
  s2All<-c(); 
  tt<-c();
  ss<-c();   
  for(j in (1:d) ){ 
    xMove=xMove1[,j] 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    if(j<=ica_k){
      ans=ans*g;
    }else{
      ans1=(ans1*sqrt((s1+s2)/n));
    }
  } 
  ret=(n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1)) -(3*n/2)*log(ans)-n*log(ans1) 
  return(-ret - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove1))) ))
}

#####################################################################
# funkcja liczaca ICA wersja I przez estymacja MLE
#####################################################################
MLE_GSGaussian_par_grad <- function(data, param, ica_k){
  grad(func=MLE_GSGaussian_par, x=param, data=data, ica_k=ica_k)
}

MLE_GSGaussian_par <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  s=matrix(param[(d+d*ica_k+1):(d+ica_k+d*ica_k)]);
  t=matrix(param[(d+ica_k+d*ica_k+1):(d+2*ica_k+d*ica_k)]);
  #for(i in (1:ica_k)){  
  #  v[,i]<-1/norm_vec(v[,i])*v[,i]
  #}  
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
  ans2 = - ans2 - n*(d/2*log(2*pi*exp(1)) + 1/2*log(det(cov(xMove))) )
  #if(ica_k == d){
  #  ans2 = - ans2 #- n*(log(1/abs(det(vv))))
  #}else{
  #  ans2 = - ans2 - n*(1/2*log(det(cov(xMove))) )    
  #}
  

  #ans2 = - ans2 - n*(log(1/abs(det(vv))))
  #ans2 = - ans2 - n*(log(det(solve((vv)%*%t(vv))%*%(vv))) ) 
  return(ans2)
}

###################################################################################################################
###################################################################################################################


################################################
# G³owna funkcja
################################################
convert_param <- function(param, x) { # param - punkt startoey # x - dane  
#  op<-optim(param, fr_inv, gr=gr_fr, x=x,method="BFGS"); # minimalizacja grafjetowa
  op<-optim(param, fr_inv, x=x); # minimalizacja grafjetowa  
  n=length(x[,1])
  d=length(x[1,])  
  c1=sqrt(2/pi)^d
  m<-matrix(op$par[1:d])
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  #v=matrix(op$par[(d+1):length(op$par)],nrow=d,byrow=TRUE)
  
#  if(det(v)<=0){
#    v<-v+diag(x = 0.000000000001, nrow=d, ncol=d)
#  }
#  vv=solve(v)
#  t<-c();
#  s<-c();
#  gAll<-c();
#  s1All<-c();
#  s2All<-c();        
  #xMove<-t(apply(data, 1, function(x) ( (solve(t(vv)%*%vv)%*%t(vv))%*%(matrix(x)-m)) ) );
#  for(j in (1:d) ){ xMove<-t(apply(x, 1, function(x) (t(vv[j,])%*%(matrix(x)-m))[1,1] )); s1<-sum(xMove[xMove<=0]^2); s1=max(s1,0.000000000001); s2<-sum(xMove[xMove>0]^2); s2=max(s2,0.000000000001); s1All<-c(s1All,s1); s2All<-c(s2All,s2); g=(s1^(1/3)+s2^(1/3)); gAll<-c(gAll,g); lTemp=( (s1)^(2/3)*g/n); lTemp=sqrt(lTemp); tTemp=((s2/s1)^(1/3)); t<-c(t,tTemp); s<-c(s,lTemp);};
#  if(op$value<=0 | is.na(op$value) | !is.finite(100)){
#    op$value=0.000000000001
#  }
#  MLE=(n*d/2)*log((2*n)/(pi*exp(1))) + (-3*n/2)*( log( op$value ));    
  
#  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=v))
}  
################################################
# funkcja, ktora minimalizujemy
################################################
fr_inv <- function(param, x) {  
  n=length(x[,1]); 
  d=length(x[1,]);  
  m<-matrix(param[1:d]); 
  
  w=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  
  w=w+diag(x = 0.000000000001, nrow=d, ncol=d)
  vv=w;
  ans=1; 
  for(j in (1:d) ){ 
    xMove<-t(apply(x, 1, function(x) (t(vv[j,])%*%(matrix(x)-m))[1,1] )); 
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3)); 
    if(is.finite( ans*g )==T){
      ans=ans*g;
    }
    if(is.finite( ans )==F){
      print(j)
      print(ans)
      print(log(ans*abs(abs(det(w))^(-2/3))))
      break;
    }
  }; 
  if(det(w)<=0){
    w<-w+diag(x = 0.000000000001, nrow=d, ncol=d)
  }
  
  ans=ans*abs(abs(det(w))^(-2/3)); 
  
  
  
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

###############################################################################################
###############################################################################################

################################################
# ICA_FEW_I
################################################
convert_param_ica_few_I <- function(x, ica_k, n_starts=1) { # param - punkt startoey # x - dane  
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  covariance<-cov(x)
  e<-eigen(covariance)
  vec_sel<-c(1,2)
  param<-c(m,e$vectors[,vec_sel])
  #op<-optim(par=param, fn=GSGaussian_par, gr=GSGaussian_par_grad, data=X, ica_k=ica_kk, method="BFGS")#, method="SANN")
  op<-optim(par=param, fn=GSGaussian_par, data=x, ica_k=ica_k,control = list(maxit = 5))
  param_temp=op$par
  value_temp=op$value
  for(i in (1:n_starts)){
    vec_sel<-sample((1:d), ica_k, replace = FALSE)
    vec_sel<-c(2,3)
    param<-c(m,e$vectors[,vec_sel])  
    print("2")
    op<-optim(par=param, fn=GSGaussian_par, data=x, ica_k=ica_k,control = list(maxit = 5))
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }
  m<-matrix(param_temp[1:d]); 
  v=matrix(param_temp[(d+1):(d+d*ica_k)],ncol=ica_k,byrow=F);
  move_data<-sweep(x, 2, m)
  return(list(m=m,v=v,par=param_temp,X=move_data))
}  

################################################
# ICA_FEW_II
################################################
convert_param_ica_few_II <- function(x, ica_k, n_starts=2) { # param - punkt startoey # x - dane 
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  covariance<-cov(x)
  e<-eigen(covariance)
  param<-c(m,e$vectors)
  #op<-optim(par=param, fn=GSGaussian_par_mle_new, gr=GSGaussian_par_mle_new_grad, data=X, ica_k=ica_kk, method="BFGS")#, method="SANN")
  op<-optim(par=param, fn=GSGaussian_par_mle_new, data=x, ica_k=ica_k,control = list(maxit = 5)) 
  param_temp=op$par
  value_temp=op$value
  for(i in (1:n_starts)){
    vec_sel<-sample((1:d), d, replace = FALSE)
    vec_sel<-c(3,2,1)
    param<-c(m,e$vectors[,vec_sel])  
    print("2")
    op<-optim(par=param, fn=GSGaussian_par_mle_new, data=x, ica_k=ica_k,control = list(maxit = 5))
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }
  n=length(x[,1]); 
  d=length(x[1,]);  
  m<-matrix(op$par[1:d]); 
  v=matrix(op$par[(d+1):(d+d*d)],ncol=d,byrow=F);
  move_data<-sweep(x, 2, m)
  
  return(list(m=m,v=v,par=op$par,X=move_data))
}  
