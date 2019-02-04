
####################step 1

###########fit Gaussina mixture model and obtain k based on BIC cretia


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
      kp<-k
    } else
    {
      k<-min(k,kmax-1)
      kp<-k+1
    }
  } else
  {
    k=kprior
    kp<-k
  }
  gmm<-densityMclust(x2,G=k,modelNames="V")
#  plot(gmm,what="density",x2,breaks=50)
  aa<-rep(0,ncol(msset))

  aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-apply(gmm$z,1, function (x) which(x==max(x)))
  msset$gmm<-aa
#  image(msset, formula = gmm~x*y,main=paste0("feature",f))
  k<-kp
  return(list(gmm,k))
}



########################
#################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

