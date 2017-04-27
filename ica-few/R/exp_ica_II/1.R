rm(list = ls())
source("F:\\AfCEC\\artykuly\\ACEC\\R\\OSTATECZNY_MODEL\\GeneralSplitGausian.R")
source("D:\\Dropbox\\ica-few\\R\\MODEL_OSTATECZNY\\ica_few_MLE_II_grad.R")

set.seed(12345)
#install.packages("heavy")
#install.packages("LambertW")
#install.packages("stochvol")
library("EMMIXuskew")
data("Lympho")

data=Lympho[sample(1:nrow(Lympho), 500), ]
#data=Lympho
x<-data[,1:2]

x<-data
m=apply(x,2,mean)
s<-rep(1,length(m))
covariance<-cov(x)
e<-eigen(covariance)

param<-c(m,t(e$vectors))

model<-convert_param(param=param, x=x) 
#model<-convert_param_gr(param=param, x=X)
model$MLE
#model$cov
#model1$cov


plot(x, main = "",pch=20,asp=1,col="azure4",xlim=c(2,6),ylim=c(2,6),xlab="",ylab="")
segments(model$m[1], model$m[2], model$m[1]+model$cov[1,1], model$m[2]+model$cov[2,1], col= 'black',lwd=3)
segments(model$m[1], model$m[2], model$m[1]+model$cov[1,2], model$m[2]+model$cov[2,2], col= 'black',lwd=3)

segments(model$m[1], model$m[2], model$m[1]-model$cov[1,1], model$m[2]-model$cov[2,1], col= 'black',lwd=3)
segments(model$m[1], model$m[2], model$m[1]-model$cov[1,2], model$m[2]-model$cov[2,2], col= 'black',lwd=3)



model_1<-convert_param_few_II(x=x, ica_k=2, n_starts=2) 

segments(model_1$m[1], model_1$m[2], model_1$m[1]+model_1$cov[1,1], model_1$m[2]+model_1$cov[2,1], col= 'red',lwd=2)
segments(model_1$m[1], model_1$m[2], model_1$m[1]+model_1$cov[1,2], model_1$m[2]+model_1$cov[2,2], col= 'red',lwd=2)

segments(model_1$m[1], model_1$m[2], model_1$m[1]-model_1$cov[1,1], model_1$m[2]-model_1$cov[2,1], col= 'red',lwd=2)
segments(model_1$m[1], model_1$m[2], model_1$m[1]-model_1$cov[1,2], model_1$m[2]-model_1$cov[2,2], col= 'red',lwd=2)


#dev.copy(pdf,'C:\\Users\\Przemek\\Desktop\\Dropbox\\ICA_gsnd\\piszemy\\images1\\distribution\\plot_2d_split_gaussian_lymp.pdf')
#dev.off()

ica_few_II(data=x,param=param,ica_k=2)
fr_inv(param=param, x=x)

gr_ica_few_II(data=x, param=param, ica_k=2)
gr_fr(param=param, x=x)

