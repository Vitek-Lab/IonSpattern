###############test model using mouse uterine data

###############calculate weight matrix

w<-W_matrix(msset,sp_ratio=3,radius=1)


ap<-rowMeans(spectra(msset))
featurelist<-which(ap>sort(ap,decreasing = T)[11])
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
    k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, initialization="km")
  }
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=2,step=1e5,Annl=0, initialization="km")
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
