
####################step 1

###########fit Gaussina mixture model and obtain k based on BIC cretia

library(Cardinal)
library("mclust")
GMM<-function(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
{
  
  lod<-0
  x<-spectra(msset)[f,]
  ###leave out pixels with intensity lower than LOD and outliers
  if (length(x[x==min(x)])/ncol(msset)>0.01)
  {
    lod<-1
    x[x==min(x)]<-NA
    x2<-x[x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2])]
    x2<-x2[!is.na(x2)]

  } else
  {
    x2<-x[x<quantile(x)[3]+out*(quantile(x)[4]-quantile(x)[2])]
  }
  #####
  if (kprior==0)
  {
    gmm<-densityMclust(x2,modelNames="V")
    if (length(gmm$BIC)>2)
    {
      Di<-rep(NA,length(gmm$BIC)-1)
      for (i in 2:length(gmm$BIC))
      {
        Di[i-1]<-gmm$BIC[i]-gmm$BIC[i-1]
      }
      if (length(which(abs(Di)<3))!=0)
      {
        k<-min(which(abs(Di)<3))
      } else
      {
        k<-gmm$G
      }
    } else
    {
      k<-gmm$G
    }
    if (lod==0)
    {
      k<-min(k,kmax)
    } else
    {
      k<-min(k,kmax-1)
      k<-k+1
    }
  } else
  {
    k=kprior
  }
  gmm<-densityMclust(x2,G=k,modelNames="V")
  plot(gmm,what="density",x2,breaks=50)
  aa<-rep(0,ncol(msset))

  aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-apply(gmm$z,1, function (x) which(x==max(x)))
  msset$gmm<-aa
  image(msset, formula = gmm~x*y,main=paste0("feature",f))
  
  return(list(gmm,k))
}



########################
#################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#################mouse brain

load("pnnl_peaksTIC-new.rdata")
pnnlCropped<-pnnl.peaksTIC_1
rm(pnnl.peaksTIC_1)
selectedSamples <- c("c1_1", "c1_2", "c1_3", 
                     "z6_1", "z6_2", "z6_3")

pnnlCropped <- pnnlCropped[,pnnlCropped$sample %in% selectedSamples]
pnnlCropped$mouse <- factor(substr(pnnlCropped$sample, 1, 2))
ap<-rowMeans(spectra(pnnlCropped))
featurelist<-which(ap>sort(ap,decreasing = T)[41])
ncomp<-matrix(0,nrow=nrow(pnnlCropped),ncol=6)
cm<-list()
diff<-rep(0,nrow(pnnlCropped))
kmax=3
out=3
for (i in 1:nrow(pnnlCropped))
{
  f<-featurelist[i]
  k<-rep(1,6)
  for (s in 1:6)
  {
    msset<-pnnlCropped[,pnnlCropped$sample==unique(pnnlCropped$sample)[s]]
    gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out)
    ncomp[f,s]<-gmm$G
  
  }
  if (length(unique(ncomp[f,]))!=1)
  {
    diff[i]=1
  }
  k=Mode(k)
  cm[[f]]<-matrix(0,ncol=k,nrow=6)
  for (s in 1:6)
  {
    msset<-pnnlCropped[,pnnlCropped$sample==unique(pnnlCropped$sample)[s]]
    gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,k_prior=k)
    cm[[f]][s,]<-gmm$parameters$mean
    
  }
  
}


