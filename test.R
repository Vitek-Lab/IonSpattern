###############test model using mouse uterine data

###############calculate weight matrix
load("data/S1RR2.rdata")
sp_ratio=1
w<-W_matrix(msset,sp_ratio=sp_ratio,radius=1)


ap<-rowMeans(spectra(msset))
featurelist<-which(ap>sort(ap,decreasing = T)[101])
print(featurelist)
mu_list<-list()
sig_list<-list()
label_list<-list()
for (i in 1:length(featurelist))
{
  f<-featurelist[i]
  k<-0
  if (k==0)
  {
    gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
    k<-gmm[[2]]
    gmm<-gmm[[1]]
    k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, step=1e5,initialization="km")
  }
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=sp_ratio,step=1e5,iteration=1000,Annl=0, sig=0.2,initialization="km")
  mu_list[[i]]<-dgmm[[1]]
  sig_list[[i]]<-dgmm[[2]]
  label_list[[i]]<-dgmm[[10]]
  
}


for (i in 1:length(featurelist))
{
  xx<-rep(1,ncol(msset))
  if (length(rmlist)!=0)
  {
    xx[-rmlist]<-apply(label_list[[i]],1, function (x) which(x==max(x))) 
  } else
  {
    xx<-apply(label_list[[i]],1, function (x) which(x==max(x)))
  }
  
  msset$dgmm<-xx
  image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F)
}



xx<-apply(label_list[[1]],1, function (x) which(x==max(x)))

xx[xx %in% c(1,2,3)]<-0

corr<-rep(0,nrow(msset))
for ( i in 1:nrow(msset))
{
  corr[i]<-cor(xx, spectra(msset)[i,])
}


pdf("corr2.pdf")
par(mfrow=c(2,2))
for ( f in 1:length(mz(msset)))
{
  image(msset,feature=f,main=paste0(f, "m/z=", round(mz(msset)[f]), "corr=", round(corr[f],digits=2)))
}
dev.off()



x<-spectra(msset)[f,]
x[x==min(x)]<-NA
x2<-x[x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2])]
x2<-x2[!is.na(x2)]
km<-kmeans(x2,centers=3)

aa<-rep(0,ncol(msset))

aa[!is.na(x) & (x<quantile(x, na.rm=TRUE)[3]+out*(quantile(x,na.rm=TRUE)[4]-quantile(x,na.rm=TRUE)[2]))]<-km$cluster
msset$dgmm<-aa
image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F)