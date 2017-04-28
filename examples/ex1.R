#rm(list = ls())
library(pixmap)
library(png)
library(ICAA)

script.dir <- dirname(sys.frame(1)$ofile)

ima1 <- readPNG(paste(script.dir, "/pics/bird.png", sep = ""))
ima2 <- readPNG(paste(script.dir, "/pics/cowboy.png", sep = ""))

# number of rows
size<-dim(ima1)[1]
p1 <- pixmapRGB(ima1, nrow=size, cellres=1)
p2 <- pixmapRGB(ima2, nrow=size, cellres=1)

data1<-c( getChannels(p1, colors = "red") )
data2<-c( getChannels(p2, colors = "red") )

# mixing images
S <- cbind(data1, data2)
A <- matrix(c(1, 1, 2, -1), 2, 2)
X <- (S %*% A)

mix1<-matrix(X[,1], nrow=size,byrow=FALSE)
imageMix1 <- pixmapRGB(mix1, nrow=size,cellres=1)
mix2<-matrix(X[,2], nrow=size,byrow=FALSE)
imageMix2 <- pixmapRGB(mix2, nrow=size,cellres=1)

plot(imageMix1)
invisible(readline(prompt="Plotting image mix 1. Press [enter] to continue.."))

plot(imageMix2)
invisible(readline(prompt="Plotting image mix 2. Press [enter] to continue.."))

# ICAA magic here -- original images are retrieved!
Snew<-ICAA( X )[[1]]

plot(pixmapRGB(matrix(Snew[,1], nrow=size,byrow=FALSE), nrow=size, cellres=1),asp=1)
invisible(readline(prompt="Plotting retrieved image 1. Press [enter] to continue.."))

plot(pixmapRGB(matrix(-Snew[,2], nrow=size,byrow=FALSE), nrow=size, cellres=1),asp=1)
invisible(readline(prompt="Plotting retrieved image 2. Press [enter] to continue.."))
