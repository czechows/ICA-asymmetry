rm(list = ls())
library(jpeg)
library(pixmap)
library(png)
library(fastICA)
library(ica)
library(PearsonICA)
library(ProDenICA)
library(steadyICA)
library(fastICA)
library(tsBSS)

components=2

#source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE.R")
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE_II_grad.R")

path<-"F:\\image_dataset\\R_1\\"
path_res<-"F:\\image_dataset\\R_reduce_ICA_images\\"
#sizes<-c(12,4,8)
pathLocal<-paste(path, 1, "\\", "", sep = "")
name=paste(1, "_", 2, sep = "")
ima1 <- readPNG(paste(pathLocal, 1, ".png", sep = ""))
ima2 <- readPNG(paste(pathLocal, 2, ".png", sep = ""))


size<-dim(ima1)[1]
p1 <- pixmapRGB(ima1, nrow=size)
p2 <- pixmapRGB(ima2, nrow=size)
canal1<-getChannels(p1, colors = "red")
canal2<-getChannels(p2, colors = "red")

data1<-c(canal1)
data2<-c(canal2)


################################################
S <- cbind(data1, data2, rnorm(n=length(data2),m=0,sd=1))#rep(1,length(data1)))
A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
#A <- matrix(c(1, 1, 0, 1, -1, 0, 0, 0, 1), 3, 3)
X <- S %*% A
################################################

#S <- cbind(data1, data2)
#A <- matrix(c(1, 1, 1, -1), 2, 2)
#X <- S %*% A
#data3 = rnorm(n=length(X[,1]),mean=0,sd=1)
#X <- cbind(X,data3)


mix1<-matrix(X[,1], nrow=size,byrow=FALSE)
imageMix1 <- pixmapRGB(mix1, nrow=size)    

mix2<-matrix(X[,2], nrow=size,byrow=FALSE)
imageMix2 <- pixmapRGB(mix2, nrow=size)    

#write.table(mix1, file = paste(path_res, name, "_", "sum",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#write.table(mix2, file = paste(path_res, name, "_", "div",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)

a11 <- icafast(X, 3, fun="logcosh")
a11$S
a11$M

###############################################

can_1=matrix(a11$S[,1],nrow=size)
can_2=matrix(a11$S[,2],nrow=size)
can_3=matrix(a11$S[,3],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1")
plot(pixmapRGB(can_2,nrow=size),main="can 2")
plot(pixmapRGB(can_3,nrow=size),main="can 3")
#image(matrix(move_data%*%((covTem))[,2],nrow=size))
#write.table(matrix(move_data%*%((covTem))[,1], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "ICA11_1",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#write.table(matrix(move_data%*%((covTem))[,2], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "ICA11_2",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


#model_1<-convert_param_few_II_3(x=X, ica_k=2, n_starts=1)
#m=apply(X,2,mean)
#A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
#param<-c(m,solve(A))
#ica_kk=2
#round(grad(func=ica_few_II, x=param, data=X, ica_k=ica_kk),digits=7)
#round(gr_ica_few_II(param=param, data=X, ica_k=ica_kk),digits=7)

model_1<-convert_param_few_II_SIPLE(x=X, ica_k=2, n_starts=1) 
convert_param_few_II_2(x=X, ica_k=2, n_starts=1) 
AA<- solve(A)
d=length(X[1,])
#m=rep(0,d)
m=apply(X,2,mean)
param1<-c(m,AA)
ica_few_II(param=param1, data=X, ica_k=2)

m<-model_1$m
v<-model_1$cov
move_data<-sweep(X, 2, m)

covT<-v
#covT[,1]<-1/norm_vec(covT[,1])*covT[,1]
#covT[,2]<-1/norm_vec(covT[,2])*covT[,2]
#covT[,3]<-1/norm_vec(covT[,3])*covT[,3]

dataSGD<-move_data%*%((covT))    

can_1=matrix(dataSGD[,1],nrow=size)
can_2=matrix(dataSGD[,2],nrow=size)
can_3=matrix(dataSGD[,3],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1 ica_few")
plot(pixmapRGB(can_2,nrow=size),main="can 2 ica_few")
plot(pixmapRGB(can_3,nrow=size),main="can 3 ica_few")

AA<-solve(A)
covT<-AA
covT[,1]<-1/norm_vec(covT[,1])*covT[,1]
covT[,2]<-1/norm_vec(covT[,2])*covT[,2]
covT[,3]<-1/norm_vec(covT[,3])*covT[,3]

dataSGD<-move_data%*%((covT))    

can_1=matrix(dataSGD[,1],nrow=size)
can_2=matrix(dataSGD[,2],nrow=size)
can_3=matrix(dataSGD[,3],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1 A")
plot(pixmapRGB(can_2,nrow=size),main="can 2 A")
plot(pixmapRGB(can_3,nrow=size),main="can 3 A")

###############################################
congru11_1_a<-congru(data1,dataSGD[,1])
congru11_1_b<-congru(data1,dataSGD[,2])
congru11_1_c<-congru(data1,dataSGD[,3])
congru11_1=max(abs(congru11_1_a),abs(congru11_1_b),abs(congru11_1_c))
congru11_1
congru11_2_a<-congru(data2,dataSGD[,1])
congru11_2_b<-congru(data2,dataSGD[,2])
congru11_2_c<-congru(data2,dataSGD[,3])
congru11_2=max(abs(congru11_2_a),abs(congru11_2_b),abs(congru11_2_c))
congru11_2
###############################################







plot3d(X)

m<-model_1$m
v1<-model_1$cov[1,1:3]
v2<-model_1$cov[2,1:3]
v3<-model_1$cov[3,1:3]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="green", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="red", pch=20, lwd=2)

pairs(X,pch=20)




v1<-a11$W[1,1:3]
v2<-a11$W[2,1:3]
v3<-a11$W[3,1:3]

segments3d(x=as.vector(t(c(m[1],m[1]+v1[1]))),
           y=as.vector(t(c(m[2],m[2]+v1[2]))),
           z=as.vector(t(c(m[3],m[3]+v1[3]))),col="blue", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v2[1]))),
           y=as.vector(t(c(m[2],m[2]+v2[2]))),
           z=as.vector(t(c(m[3],m[3]+v2[3]))),col="blue", pch=20, lwd=4)

segments3d(x=as.vector(t(c(m[1],m[1]+v3[1]))),
           y=as.vector(t(c(m[2],m[2]+v3[2]))),
           z=as.vector(t(c(m[3],m[3]+v3[3]))),col="blue", pch=20, lwd=2)



#model_X<-move_data
#for(i in (1:length(covTem[1,]))){
#  covTem[,i]<-1/norm_vec(covTem[,i])*covTem[,i]
#}  
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,1],nrow=size),nrow=size),main=file_name)
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,2],nrow=size),nrow=size),main=file_name)    
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,3],nrow=size),nrow=size),main=file_name)