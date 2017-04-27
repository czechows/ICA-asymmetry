#install.packages('pixmap')
#install.packages('png')
#install.packages('jpeg')
#install.packageslibrary('fastICA')
#install.packages('seewave')

source("F:\\AfCEC\\artykuly\\ACEC\\R\\OSTATECZNY_MODEL\\GeneralSplitGausian.R")
source("F:\\AfCEC\\artykuly\\ACEC\\R\\OSTATECZNY_MODEL\\LogisticDistribution.R")

library("fastICA")
library("ica")

set.seed(1222)

rlogistic<- function(x, m, s,size){
  n <- 400000
  X <- runif(n,min=-30,max=30)  # generowanie z rozk³adu jednostajnego na [-1,1]
  Y <- runif(n,min=0,max=1)
  ans<-c()
  for(i in (1:length(X)) ){
    if( Y[i] < logistic1D(x=X[i], m=m, s=s)  ){
      ans<-c(ans,X[i])
    }
  }
  return(ans[1:size])
}

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


doICA<-function(vawe1_l, vawe2_l, j){
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
  param<-c(m,e$vectors)
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
  
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_SG_a_",j,".pdf",sep=""))
  #dev.off()
  
  plot(X, main = "SG",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  #  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  #  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  #  dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_SG_b_",j,".pdf",sep=""))
  #  dev.off()
  
  plot(X, main = "SG",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  #  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_SG_c_",j,".pdf",sep=""))
  #dev.off()
  
  plot(X, main = "SG",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  #  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  #  dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_SG_d_",j,".pdf",sep=""))
  #  dev.off()
  
  
  #play(normalize(w_dataSGD_1,unit = "8"))
  #play(normalize(w_dataSGD_2,unit = "8"))
  
  
  ###############################################
  ###############################################
  a11 <- icafast(X, 2, fun="logcosh")
  a11$S
  a11$M
  
  m=c(0,0)
  model_cov=a11$M
  model_cov[,1]<-1/norm_vec(model_cov[,1])*model_cov[,1]
  model_cov[,2]<-1/norm_vec(model_cov[,2])*model_cov[,2]
  
  plot(X, main = "FastICA",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  #  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_ica_11_a_",j,".pdf",sep=""))
  #dev.off()
  
  plot(X, main = "FastICA",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  #segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  #segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_ica_11_b_",j,".pdf",sep=""))
  #dev.off()
  
  plot(X, main = "FastICA",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  #segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  #segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_ica_11_C_",j,".pdf",sep=""))
  #dev.off()
  
  plot(X, main = "FastICA",pch=20,asp=1,col="azure4",xlim=c(-3.5,3.5),ylim=c(-3.5,3.5),xlab="",ylab="")
  segments(model$m[1], model$m[2], model$m[1]+model_cov[1,1], model$m[2]+model_cov[2,1], col= 'black',lwd=3)
  #segments(model$m[1], model$m[2], model$m[1]+model_cov[1,2], model$m[2]+model_cov[2,2], col= 'black',lwd=3)
  
  #segments(model$m[1], model$m[2], model$m[1]-model_cov[1,1], model$m[2]-model_cov[2,1], col= 'black',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]-model_cov[1,2], model$m[2]-model_cov[2,2], col= 'black',lwd=3)
  
  segments(model$m[1], model$m[2], model$m[1]-A[1,1], model$m[2]-A[2,1], col= 'red',lwd=3)
  segments(model$m[1], model$m[2], model$m[1]+A[1,2], model$m[2]+A[2,2], col= 'red',lwd=3)
  
  
  #dev.copy(pdf,paste("C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\syntetic\\plot_ica_11_D_",j,".pdf",sep=""))
  #dev.off()
  
  
  
  a12 <- icafast(X, 2, fun="exp")
  a12$S
  a12$M
  
  a13 <- icafast(X, 2, fun="kur")
  a13$S
  a13$M
  
  a21 <- icaimax(X, 2, fun="tanh" )
  a21$S
  a21$M
  
  a22 <- icaimax(X, 2, fun="log" )
  a22$S
  a22$M
  
  a23 <- icaimax(X, 2, fun="ext" )
  a23$S
  a23$M
  
  a3 <- icajade(X, 2)
  a3$S
  a3$M
  ###############################################
  ###############################################
  
  
  congru0_1<-congru(data1,dataSGD[,1])
  congru0_2<-congru(data2,dataSGD[,2])
  
  congru11_1<-congru(data1,a11$S[,1])
  congru11_2<-congru(data2,a11$S[,2])
  congru12_1<-congru(data1,a12$S[,1])
  congru12_2<-congru(data2,a12$S[,2])
  congru13_1<-congru(data1,a13$S[,1])
  congru13_2<-congru(data2,a13$S[,2])
  
  congru21_1<-congru(data1,a21$S[,1])
  congru21_2<-congru(data2,a21$S[,2])
  congru22_1<-congru(data1,a22$S[,1])
  congru22_2<-congru(data2,a22$S[,2])
  congru23_1<-congru(data1,a23$S[,1])
  congru23_2<-congru(data2,a23$S[,2])
  
  congru3_1<-congru(data1,a3$S[,1])
  congru3_2<-congru(data2,a3$S[,2])
  
  congruCompare=c(congru0_1,congru0_2,
                  congru11_1, congru11_2, congru12_1, congru12_2, congru13_1,congru13_2,
                  congru21_1, congru21_2, congru22_1, congru22_2, congru23_1,congru23_2,
                  congru3_1, congru3_2)
  
  
  acy0<-acy(A,t(solve(covT)))
  #acy0<-acy(A,covT)  
  
  acy11<-acy(A,a11$M)
  acy12<-acy(A,a12$M)
  acy13<-acy(A,a13$M)
  
  acy21<-acy(A,a21$M)
  acy22<-acy(A,a22$M)
  acy23<-acy(A,a23$M)
  
  acy3<-acy(A,a3$M)
  
  acyCompare=c(acy0,
               acy11, acy12, acy13,
               acy21, acy22, acy23,
               acy3)
  
  return(list(sgdX=dataSGD, cov=covT, stat2=congruCompare, stat3=acyCompare))
}

library(seewave)
library(tuneR)
library(fastICA)

#source("C:\\Users\\Przemek\\Desktop\\Dropbox\\ACEC\\R\\OSTATECZNY_MODEL\\GeneralSplitGausian.R")
#source("F:\\image_dataset\\R\\imageICA.R")


path<-"F:\\image_dataset\\Sytntetic_data\\"
#sizes<-c(12,4,8)


for(j in c(2)){
  if(j==1){
    m=c(0,0)
    s=c(4,3)
    vawe1_l = rlogistic(m=0,s=s[1],size=1000)
    vawe2_l = rlogistic(m=0,s=s[2],size=1000)    
  }else{
    m=c(0,0)
    s=c(1,1)
    t=c(0.5,0.2)
    vawe1_l = rspliNorm(m=0,s=s[1],t=t[1],size=1000)
    vawe2_l = rspliNorm(m=0,s=s[2],t=t[2],size=1000)
  }
  ICA_vs_SGD<-doICA(vawe1_l, vawe2_l,j)
  #write(ICA_vs_SGD$stat2,paste(path,"statistic2.txt", sep = ""), ncolumns=length(ICA_vs_SGD$stat2),append = TRUE)
  #write(ICA_vs_SGD$stat3,paste(path,"statistic3.txt", sep = ""), ncolumns=length(ICA_vs_SGD$stat3),append = TRUE)    
  
}

#write.table(data1, file = paste(path, "1.txt", sep = ""), append = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
