###########example of DGMM component wise statistical analysis: we use mouse brain data, which has 18 samples from two groups
library(nlme)
library("Biobase")
library(Cardinal)
library(spam)
library(mvtnorm)
library("mclust")
library(gridExtra)
library("ggplot2")
library(RColorBrewer)
############load data
load("data/pnnl_peaksTIC-new.rdata")
pnnlCropped<-pnnl.peaksTIC_1

##choose feature 321 of mouse brain data

f<-321

##Estimate how many Gaussian components it has for feature 321
k_list<-rep(1,length(unique(pnnlCropped$sample)))

kmax=5
out=10

for (i in 1:length(unique(pnnlCropped$sample)))
{
  
  msset<-pnnlCropped[,pnnlCropped$sample==unique(pnnlCropped$sample)[i]]
  gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
  k<-gmm[[2]]
  gmm<-gmm[[1]]
  k_list[i]<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, step=1e5,initialization="gmm",r_max=3)
  
}
k<-Mode(k_list)
######print number of components
print(paste0("number of components is: ", k))

#####initializing a matrix for storing mean intensities of each component for each sample
means<-matrix(0,ncol=k,nrow=length(unique(pnnlCropped$sample)))
for (i in 1:length(unique(pnnlCropped$sample)))
{
  msset<-pnnlCropped[,pnnlCropped$sample==unique(pnnlCropped$sample)[i]]
  sp_ratio=5
  w<-W_matrix(msset,sp_ratio=sp_ratio,radius=1)
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=sp_ratio,step=1e5,iteration=1000,Annl=0, sig=0.2,initialization="km")
  means[i,]<-dgmm[[1]]
}

###################the mean intensities are ranked in ascending order for all samples, thus in this case, components were 
###################matched among samples by their relative orders
########t-test
########test mean intensity directly
for (i in 1:k)
{
  ttest_p<-t.test(means[1:9,i],means[10:18,i])$p.value
  print(paste0("component ", i, ": ", ttest_p))
}

########test component versus component ratio

for ( i in 2:k)
{
  xx<-means[,i]/means[,1]
  ttest_p<-t.test(xx[1:9],xx[10:18])$p.value
  print(paste0("component ", i, "/component 1",  ": ", ttest_p))
}


#########linear mixed model
########test mean intensity directly

xx1<-c(rep(1:3,4),1,1:2,1:3)
xx2<-c(rep("sal",9),rep("cpg",9))
xx3<-as.character(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),5,rep(6,2),rep(7,3)))

for (i in 1:k)
{
  xx<-means[,i]
  m1 <- lme(xx~xx2,random=~1|xx3)
  print(paste0("component ", i, ": "))
  print(anova(m1))
}

########test component versus component ratio
for (i in 2:k)
{
  xx<-means[,i]/means[,1]
  m1 <- lme(xx~xx2,random=~1|xx3)
  print(paste0("component ", i, "/component 1",  ": "))
  print(anova(m1))
}

