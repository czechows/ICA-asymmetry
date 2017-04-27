library(rgl)
data=cbind(rnorm(100,0,1),rnorm(100,0,1),rnorm(100,0,1))
plot3d(data)

m=apply(data,2,mean)
covariance<-cov(data)
e<-eigen(covariance)
v1=e$vectors[,1]
v2=e$vectors[,2]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="red", pch=20)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="red", pch=20)

vv=e$vectors
x=data
j=1



xMove<-t(apply(x, 1, function(x) ( solve((vv[1:2,])%*%t(vv[1:2,]))%*%(vv[1:2,]%*%(matrix(x)-m)) ) ));
xMove=xMove[,1]
dim(xMove)
plot(xMove,pch=1)
xx<-vv[1:2,]
xMove1<-t(apply(xx, 1, function(xx) ( solve((vv[1:2,])%*%t(vv[1:2,]))%*%(vv[1:2,]%*%(matrix(xx))) ) ));
dim(xMove1)
xMove1<-rbind(xMove1,c(0,0))
points(xMove1,pch=20,col="red")

xMove<-t(apply(x, 1, function(x) (vv[j,]%*%(matrix(x)-m))[1,1] ));

points3d(t(xMove)*vv[j,],col="green", size=3)
