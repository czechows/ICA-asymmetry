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
#A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
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
m=apply(X,2,mean)
move_data<-sweep(X, 2, m)



can_1=matrix(a11$S[,1],nrow=size)
can_2=matrix(a11$S[,2],nrow=size)
can_3=matrix(a11$S[,3],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1")
plot(pixmapRGB(can_2,nrow=size),main="can 2")
plot(pixmapRGB(can_3,nrow=size),main="can 2")
#image(matrix(move_data%*%((covTem))[,2],nrow=size))
#write.table(matrix(move_data%*%((covTem))[,1], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "ICA11_1",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#write.table(matrix(move_data%*%((covTem))[,2], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "ICA11_2",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)


model_1<-convert_param_few_II(x=X, ica_k=3, n_starts=1) 

op$par
n=length(X[,1]); 
d=length(X[1,]);  
m<-matrix(op$par[1:d]); 
v=matrix(op$par[(d+1):(d+d*ica_kk)],ncol=ica_kk,byrow=F);
s=matrix(op$par[(d+d*ica_kk+1):(d+ica_kk+d*ica_kk)]);
t=matrix(op$par[(d+ica_kk+d*ica_kk+1):(d+2*ica_kk+d*ica_kk)]);

#covariance<-cov(X)
#e<-eigen(covariance)
#param<-c(m,t(e$vectors))
#model<-convert_param(param=param, x=X)  

covT<-v
covT[,1]<-1/norm_vec(covT[,1])*covT[,1]
covT[,2]<-1/norm_vec(covT[,2])*covT[,2]

dataSGD<-move_data%*%((covT))    

can_1=matrix(dataSGD[,1],nrow=size)
can_2=matrix(dataSGD[,2],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1 ica_few")
plot(pixmapRGB(can_2,nrow=size),main="can 2 ica_few")

#write.table(matrix(dataSGD[,1], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "SGD_1",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
#write.table(matrix(dataSGD[,2], nrow=size,byrow=FALSE), file = paste(path_res, name, "_", "SGD_2",".txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)  

move_data<-sweep(X, 2, m)

#####

modal_2=convert_param_few_II(x=X, ica_k=3, n_starts=1)
m<-modal_2$m
v<-modal_2$cov

covT<-v
covT[,1]<-1/norm_vec(covT[,1])*covT[,1]
covT[,2]<-1/norm_vec(covT[,2])*covT[,2]

dataSGD<-move_data%*%((covT))    

can_1=matrix(dataSGD[,1],nrow=size)
can_2=matrix(dataSGD[,2],nrow=size)
can_3=matrix(dataSGD[,3],nrow=size)
plot(pixmapRGB(can_1,nrow=size),main="can 1 ica_few")
plot(pixmapRGB(can_2,nrow=size),main="can 2 ica_few")
plot(pixmapRGB(can_3,nrow=size),main="can 3 ica_few")

#model_X<-move_data
#for(i in (1:length(covTem[1,]))){
#  covTem[,i]<-1/norm_vec(covTem[,i])*covTem[,i]
#}  
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,1],nrow=size),nrow=size),main=file_name)
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,2],nrow=size),nrow=size),main=file_name)    
#plot(pixmapRGB(matrix(model_X%*%((covTem))[,3],nrow=size),nrow=size),main=file_name)