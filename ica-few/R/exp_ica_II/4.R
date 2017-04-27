rm(list = ls())
library(fastICA)
library(ica)
library(PearsonICA)
library(ProDenICA)
library(steadyICA)
library(fastICA)
library(tsBSS)
library(seewave)
library(tuneR)

components=2
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE_II_grad.R")

path<-"F:\\image_dataset\\R_1_sound_few_I\\"
pathLocal<-paste(path, "vawe\\source", "", sep = "")
vawe1 <- readWave(paste(pathLocal, 1, ".wav", sep = ""))
vawe1_l = vawe1@left
vawe2 <- readWave(paste(pathLocal, 1+1, ".wav", sep = ""))
vawe2_l = vawe2@left

sdd=max(sd(vawe1_l),sd(vawe2_l))
S <- cbind(vawe1_l, vawe2_l, rnorm(n=length(vawe2_l),m=0,sd=sdd+2))
A <- matrix(c(1, 1, -1, 1, -1, 1, 1, -1, -1), 3, 3)
#  A <- matrix(c(1/2, 1/2, 1/2, -1/2), 2, 2)
X <- S %*% A
#plot(X)
#m=apply(X,2,mean)
#X<-sweep(X, 2, m)

vawe_sum = X[,1]
vawe_div = X[,2]

data1<-vawe1_l
data2<-vawe2_l

w_sum = Wave(vawe_sum, samp.rate = 8000, bit=8)
w_div = Wave(vawe_div, samp.rate = 8000, bit=8)
#savewav(w_sum, filename = paste(path, "vawe_save\\", j, "_", j+1, "_sum", ".wav", sep = ""))
#savewav(w_div, filename = paste(path, "vawe_save\\", j, "_", j+1, "_div", ".wav", sep = ""))

#play(normalize(w_sum,unit = "8"))
#play(normalize(w_div, unit = "8"))

m=apply(X,2,mean)
covariance<-cov(X)
e<-eigen(covariance)
param<-c(m,t(e$vectors))
print("ica_2")
model_2=convert_param_few_II_SIPLE(x=X, ica_k=2, n_starts=5)
print("ica_2")
m<-model_2$m
v<-model_2$cov
model_X<-sweep(X, 2, m)
dataSGD<-X%*%((covT))    


###############################################
congru11_1_a<-congru(data1,a11$S[,1])
congru11_1_b<-congru(data1,a11$S[,2])
congru11_1=max(abs(congru11_1_a),abs(congru11_1_b))
congru11_1
congru11_2_a<-congru(data2,a11$S[,1])
congru11_2_b<-congru(data2,a11$S[,2])
congru11_2=max(abs(congru11_2_a),abs(congru11_2_b))
congru11_2
###############################################


a11 <- icafast(X, 3, fun="logcosh")
a11$W%*%a11$M
a11$S
a11$M

###############################################
congru11_1_a<-congru((data1),a11$S[,1])
congru11_1_b<-congru(data1,a11$S[,2])
congru11_1=max(abs(congru11_1_a),abs(congru11_1_b))
congru11_1
congru11_2_a<-congru(data2,a11$S[,1])
congru11_2_b<-congru(data2,a11$S[,2])
congru11_2=max(abs(congru11_2_a),abs(congru11_2_b))
congru11_2
###############################################

AA<-solve(A)
covT<-AA
d=length(X[1,]);
m=rep(0,d)

dataSGD<-X%*%((covT))    

###############################################
congru11_1_a<-congru(data1,dataSGD[,1])
congru11_1_b<-congru(data1,dataSGD[,2])
congru11_1=max(abs(congru11_1_a),abs(congru11_1_b))
congru11_1
congru11_2_a<-congru(data2,dataSGD[,1])
congru11_2_b<-congru(data2,dataSGD[,2])
congru11_2=max(abs(congru11_2_a),abs(congru11_2_b))
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