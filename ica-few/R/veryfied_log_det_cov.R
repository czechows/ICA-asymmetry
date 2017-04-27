library("numDeriv")
norm_vec <- function(x) sqrt(sum(x^2))
len=500
x=cbind(rnorm(len,0,2),rnorm(len,0,1),rnorm(len,0,1))
d=2
m=apply(x,2,mean)
w=eigen(cov(x))$vectors
par=c(m,w)


fr <- function(param, data) {  
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  v=matrix(param[(d+1):length(param)],ncol=d,byrow=F)
  xMove1<-t(apply(data, 1, function(x) ( (solve(t(v)%*%v)%*%t(v))%*%(matrix(x)-m)) ) );
  ans = cov(xMove1)
  return(1/2*log( det(ans) ));
} 

gr_fr <- function(param, data) {  
  n=length(data[,1]); 
  d=length(data[1,]);  
  m<-matrix(param[1:d]); 
  ans1=rep(0,d);
  ans2 = -2*(solve(t(v)%*%v)%*%v) + solve( cov(data%*%t(v)) )%*%cov(data)%*%v
  #print(ans2)
  ans=c(ans1,c(ans2))
  return(ans);
} 

fr(par,x)

round(grad(func=fr, x=par, data=x),digits=7)
round(gr_fr(param=par, data=x),digits=7)


