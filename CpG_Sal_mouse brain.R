library("Biobase")
library(Cardinal)
library(spam)
library(mvtnorm)
library("mclust")
library(gridExtra)
library("ggplot2")
library(RColorBrewer)
load("data/pnnl_peaksTIC-new.rdata")
pnnlCropped<-pnnl.peaksTIC_1


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
msset<-pnnlCropped[,pnnlCropped$sample=="c1_1"]
set.seed(1234)
ssc <- spatialShrunkenCentroids(msset, r=1, k=4, s=3, method="adaptive",main="spatial centroid segmentation")
image(ssc, col=cbPalette[1:4], key=FALSE,asp=5,ylim=c(-10,22),xlim=c(-10,180))
skm <- spatialKMeans(msset, r=1, k=4, method="adaptive")
image(skm, col=cbPalette[1:4], key=FALSE,asp=5,ylim=c(-10,22),xlim=c(-10,180),main="spatial k-means")


sp_ratio=5
radius=1
w<-W_matrix(msset,sp_ratio=sp_ratio,radius=1)
f<-327
kmax=5
out=10
kprior=0
gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
k<-gmm[[2]]
gmm<-gmm[[1]]
k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, step=1e5,initialization="gmm",r_max=3)
dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=sp_ratio,step=1e5,iteration=1000,Annl=1, sig=0.2,initialization="km")

###########plot
xx<-apply(dgmm[[10]],1, function (x) which(x==max(x)))
skm@resultData$`r = 1, k = 4`$cluster<-as.factor(xx)
image(skm, col=c(cbPalette[1:4]), key=FALSE,asp=5,ylim=c(-10,22),xlim=c(-10,180),main="spatial DGMM: m/z=844.67") 

  
  
  f<-291
  kmax=5
  out=10
  kprior=0
  gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
  k<-gmm[[2]]
  gmm<-gmm[[1]]
  k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, step=1e5,initialization="gmm",r_max=3)
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=sp_ratio,step=1e5,iteration=1000,Annl=1, sig=0.2,initialization="gmm")
  
  ###########plot
  xx<-apply(dgmm[[10]],1, function (x) which(x==max(x)))
  skm@resultData$`r = 1, k = 4`$cluster<-as.factor(xx)
  image(skm, col=c(cbPalette[1:4]), key=FALSE,asp=5,ylim=c(-10,22),xlim=c(-10,180),main="spatial DGMM: m/z=756.51") 
  
  f<-321
  kmax=5
  out=10
  kprior=0
  gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
  k<-gmm[[2]]
  gmm<-gmm[[1]]
  k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, step=1e5,initialization="gmm",r_max=3)
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=sp_ratio,step=1e5,iteration=1000,Annl=1, sig=0.2,initialization="km")
  
  ###########plot
  xx<-apply(dgmm[[10]],1, function (x) which(x==max(x)))
  skm@resultData$`r = 1, k = 4`$cluster<-as.factor(xx)
  image(skm, col=c(cbPalette[1:4]), key=FALSE,asp=5,ylim=c(-10,22),xlim=c(-10,180),main="spatial DGMM: m/z=826.76") 