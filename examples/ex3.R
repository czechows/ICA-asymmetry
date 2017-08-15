#rm(list = ls())
library(ICAA)


X=read.table("F:\\image_dataset\\time\\CZE\\ICA_1_10.txt")
X
Snew<-ICAA( as.matrix(X), generalized=1, minimum = -10., accuracy = 1e-5 )
