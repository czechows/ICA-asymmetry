
norm_vec <- function(x) sqrt(sum(x^2))
library(matrixcalc)

convert_param_few_II <- function(x, ica_k, n_starts) { # param - punkt startoey # x - dane  
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  covariance<-cov(x)
  e<-eigen(covariance)
  param<-c(m,t(e$vectors))
  #A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
  #AA<- solve(A)
  #m=rep(0,d)
  #param<-c(m,(A))
  op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
  param_temp=op$par
  value_temp=op$value

  for(i in (1:1)){
    vec_sel<-sample((1:d), d, replace = FALSE)
    vec_sel<-c(3,2,1)
    param<-c(m,t(e$vectors)[,vec_sel]) 
    #param<-c(m,(solve(A[,vec_sel])))  
    print("2")
    op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }

  for(i in (1:n_starts)){
    vec_sel<-sample((1:d), d, replace = FALSE)
    #vec_sel<-c(3,2,1)   
    exampledf<- matrix(ceiling(runif(9,0,5)), ncol=3)
    while( is.singular.matrix(exampledf) == TRUE){
      exampledf<- matrix(ceiling(runif(9,0,5)), ncol=3)
    }
#    m=rep(0,d)
    param<-c(m,t(exampledf))  
    #param<-c(m,(solve(A[,vec_sel])))  
    print("3")
    op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }
  
  m<-matrix(param_temp[1:d])
  v=matrix(param_temp[(d+1):length(param_temp)],ncol=d,byrow=F);
  if(det(v)<=0){
    v<-v+diag(x = 0.1e-320, nrow=d, ncol=d)
  }
  vv=solve(v)    
  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=t(v),par=param_temp,X=move_data))
}  

convert_param_few_II_2 <- function(x, ica_k, n_starts, pca=0) { # param - punkt startoey # x - dane  
  if(pca>0){
    fit <- princomp(x, scores=T, cor=TRUE)
    x=fit$scores[,1:pca]
    ica_k<-min(ica_k,ica_k)
  }
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  covariance<-cov(x)
  e<-eigen(covariance)
  param<-c(m,t(e$vectors)) 
  op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS");#,control = list(maxit = 40) # minimalizacja grafjetowa
  param_temp=op$par
  value_temp=op$value
  
  m<-matrix(param_temp[1:d])
  v=matrix(param_temp[(d+1):length(param_temp)],ncol=d,byrow=F);
  #if(det(v)<=0){
  #  v<-v+diag(x = 0.1e-320, nrow=d, ncol=d)
  #}
  #vv=solve(v)    
  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=v,par=param_temp,X=move_data))
}  

convert_param_few_II_3 <- function(x, ica_k, n_starts) { # param - punkt startoey # x - dane  
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  covariance<-cov(x)
  e<-eigen(covariance)
  param<-c(m,t(e$vectors))
  #A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
  #AA<- solve(A)
  #m=rep(0,d)
  #param<-c(m,AA)
  op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
  param_temp=op$par
  value_temp=op$value
  
  for(i in (1:10)){
    vec_sel<-sample((1:d), d, replace = FALSE)
    #vec_sel<-c(3,2,1)
    m=apply(x,2,mean)
    param<-c(m,t(e$vectors[,vec_sel])) 
    #param<-c(m,(solve(A[,vec_sel])))  
    print("2")
    op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }

  for(i in (1:n_starts)){
    vec_sel<-sample((1:d), d, replace = FALSE)
    #vec_sel<-c(3,2,1)   
    exampledf<- matrix(ceiling(runif(9,0,5)), ncol=3)
    while( is.singular.matrix(exampledf) == TRUE){
      exampledf<- matrix(ceiling(runif(9,0,5)), ncol=3)
    }
    param<-c(m,t(exampledf))  
    print("3")
    op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
    if( op$value < value_temp){
      param_temp=op$par
      value_temp=op$value
    }
  }
  
  
  m<-matrix(param_temp[1:d])
  v=matrix(param_temp[(d+1):length(param_temp)],ncol=d,byrow=F);
  if(det(v)<=0){
    v<-v+diag(x = 0.1e-320, nrow=d, ncol=d)
  }
  vv=solve(v)    
  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=(v),par=param_temp,X=move_data))
}  

convert_param_few_II_SIPLE <- function(x, ica_k, n_starts) { # param - punkt startoey # x - dane  
  n=length(x[,1]); 
  d=length(x[1,]);
  m=apply(x,2,mean)
  library(fastICA)
  #a11 <- icafast(X, 3, fun="logcosh")
  #vec_sel<-c(3,2,1)
  #AA=a11$W
  A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
  AA<- solve(A)
  #m=rep(0,d)
  param<-c(m,AA) 
  op<-optim(par=param, fn=ica_few_II, gr=gr_ica_few_II, data=x, ica_k=ica_k, method="BFGS"); # minimalizacja grafjetowa
  param_temp=op$par
  value_temp=op$value
  
  m<-matrix(param_temp[1:d])
  v=matrix(param_temp[(d+1):length(param_temp)],ncol=d,byrow=F);
  move_data<-sweep(x, 2, m)
  return(list(m=m,cov=v,par=param_temp,X=move_data))
}  


ica_few_II <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  ans=1
  ans1=1
  vv=v
  #xMove1<-t(apply(data, 1, function(x) ( t(vv)%*%(matrix(x)-m)) ) );  
  for(j in (1:d) ){ 
    #xMove=xMove1[,j] 
    xMove<-t(apply(data, 1, function(x) (vv[,j]%*%(matrix(x)-m))[1,1] ));
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    g=(s1^(1/3)+s2^(1/3));
    if(j<=ica_k){
      ans=ans*g;
    }else{
      ans1=(ans1*(s1+s2)/n);
      #print("...")
      #print(ans1)
      #print(s1+s2)
      #print("...")
    }    
  } 
  #param_temp<-c(m,v,ss,tt)
  #cc<-(n*(q-p/2))*log(2)+(q*n/2)*log(n)-p*n/2*log(pi*exp(1))
  #ret=cc+(-3*n/2)*(log(ans1^(1/3))+log(ans)+log(abs(det(v))^(2/3)))
  ret=( (ans1^(1/3))*(ans)*(abs(det(vv))^(-2/3))  )
  #ret=(log(ans1^(1/3))+log(ans)+log(abs(det(v))^(2/3)))
  print("...1")
  #print(ans)
  #print(ans1)
  print(ret)
  print(log(ret))
  print("...2")
  if(ret<=0){
    return(log(0.1e-320));
  }else{
    if(is.infinite(log(ret))==T | is.na(log(ret))==T){
      print("Something went wrong")
      print("The valu of cost functio is")
      print(ret)
      print("The dimention of data is to large")
      break
    }
    return(log(ret))
  }
}




gr_ica_few_II <- function(data, param, ica_k){
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):(d+d*d)],ncol=d,byrow=F);
  #v=matrix(param[(d+1):length(param)],nrow=d,byrow=TRUE)
  #vv=(solve(t(v)%*%v)%*%t(v))
  vv=v
#  xMove1<-t(apply(data, 1, function(x) (t(vv)%*%(matrix(x)-m))));
  ans1=rep(0,d);   
  ans_c<-c()
  for(j in (1:d) ){ 
#    xMove=xMove1[,j] 
    xMove<-t(apply(data, 1, function(x) (vv[,j]%*%(matrix(x)-m))[1,1] ));
    xMove_N<-t(apply(data, 1, function(x) ((vv[,j]%*%(matrix(x)-m))[1,1] *vv[,j]) ));     
    xMove_N_w<-t(apply(data, 1, function(x) (2*(vv[,j]%*%(matrix(x)-m))[1,1]*(matrix(x)-m) ) ));            
##### old    
#    xMove_N<-t(apply(data, 1, function(x) ( (t(vv)%*%(matrix(x)-m))[j,]*t(vv)[j,] ) ) );     
#    xMove_N_w<-t(apply(data, 1, function(x) ( (2*t(vv)%*%(matrix(x)-m))[j,]*(matrix(x)-m) ) ) );  
##### old    
    
    s1<-sum(xMove[xMove<=0]^2); 
    s2<-sum(xMove[xMove>0]^2); 
    s1_N<-apply(xMove_N[xMove<=0,], 2, sum) 
    s2_N<-apply(xMove_N[xMove>0,], 2, sum) 
    s1_N_w<-apply(xMove_N_w[xMove<=0,], 2, sum) 
    s2_N_w<-apply(xMove_N_w[xMove>0,], 2, sum)
    
    if(j<=ica_k){
      g=(s1^(1/3)+s2^(1/3));     
      gg=(2/3)*(s1^(-2/3))*s1_N+(2/3)*(s2^(-2/3))*s2_N
      ggg=(1/3)*(s1^(-2/3))*s1_N_w+(1/3)*(s2^(-2/3))*s2_N_w
      ans1=ans1-1/g*gg;
      ans_c<-c(ans_c,1/g*ggg)
    }else{
      gg_1=(2/(3*(s1+s2)))*(s1_N+s2_N)
      ggg_1=(1/(3*(s1+s2)))*(s1_N_w+s2_N_w)       
      ans1=ans1-gg_1;
      ans_c<-c(ans_c,ggg_1)
    }    
  }
  if(det(vv)<=0){
    vv<-vv+diag(x = 0.1e-320, nrow=d, ncol=d)
  }
  ans=c(ans1, ans_c+ c((-2/3)*(solve(t(vv)))) )
  print("...3")
  print(ans)
  print("...4")
  return(ans);
}

