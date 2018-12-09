
#################################fit Dirichlet Gaussian Mixture Model


fit_DGMM<-function(w=w,msset=msset, f=f, k=0,out=3, kmax=3)
{
  ########## k is not specified
  if (k==0)
  {
    gmm<-GMM(msset=msset,f=f,kmax=kmax,out=out,kprior=0)
    k<-gmm[[2]]
    gmm<-gmm[[1]]
    k<-K_DGMM(msset=msset,gmm=gmm,f=f,k=k, initialization="km")
  }
  dgmm<-DGMM(msset=msset, w<-w, k=k, f=f, sp_ratio=2,step=1e5,Annl=0, initialization="Km")
  return(dgmm)
}


#############class comparison: select differentially abundant ions

############ compare = c("t-test", "ROC")
ClassCompare<-function(msset=msset, compare="t-test", f=f, k=0,out=3, kmax=3)
{
 
  if (compare=="t-test")
  {
    ###############average t-test
    meansNormal <- sapply(levels(factor(msset$sample)),function(s){
      rowMeans(spectra(msset[,msset$sample == s& msset$diagnosis== unique(msset$diagnosis)[1] ]))
    }
    )
    
    meansCancer <- sapply(levels(factor(msset$sample)),function(s){
      rowMeans(spectra(msset[,msset$sample == s& msset$diagnosis == unique(msset$diagnosis)[2]]))
    }
    )
    
    
    fitsUnpaired <- lapply(1:nrow(msset), function(i){t.test(meansCancer[i,], meansNormal[i,], paired = F)})
    pvalueUnpaired <- unlist(lapply(fitsUnpaired, function(x) x$p.value))
    av_p<-p.adjust(pvalueUnpaired, method = "BH") 
    # table(pvalueUnpaired < .05)
    
    ##############per pixel p value
    fitsUnpaired <- lapply(1:nrow(msset), function(i){t.test(spectra(msset[i,msset$diagnosis==unique(msset$diagnosis)[1]]), spectra(msset[i,msset$diagnosis==unique(msset$diagnosis)[1]]), paired = F)})
    pvalueUnpaired <- unlist(lapply(fitsUnpaired, function(x) x$p.value))
    pixel_p<-p.adjust(pvalueUnpaired, method = "BH") 
    #table(pvalueUnpaired < .05)
    
    ##############sub-tissue morphological t-test
    est_m<-list()
    for (f in 1:length(mz(msset)))
    {
      est_m[[f]]<-matrix(0,ncol=kmax,nrow=length(unique(msset$sample)))
      for ( i in 1:length(unique(msset$sample)))
      {
        dgmm<-fit_DGMM(msset=msset[,msset$sample=unique(msset$sample)[i]], f=f, k=0,out=3, kmax=3)
        est_m[[f]][i,]<-dgmm$mu
      }
    }
    
    
    x=c(1,2,3)
    pp<-matrix(0,ncol=kmax,nrow=40)
    for (f in 1:length(mz(msset)))
    {
      k<-length(x[colSums(est_m[[f]])!=0])
      for (i in 1:k)
      {
        pp[f,i]<-t.test(est_m[[f]][c(1,3,5,7,9),i], est_m[[f]][c(2,4,6,8,10),i])$p.value
      }
    }
    pp[pp==0]=NA
    pp<-matrix(p.adjust(pp, method = "BH"),ncol=kmax,byrow=FALSE)
    dgmm_p<-pp
    
    
    
    ###########plot
    
    data<-data.frame(av_p,pixel_p,dgmm_p)
    colnames(data)<-c("Average", "Per_Pixel", 'DGMM_comp1', "DGMM_comp2", "DGMM_comp3")
    data$mz<-1:40
    group<-c(rep("Averge",40),rep("Per_Pixel",40),rep("DGMM_comp1",40),rep("DGMM_comp2",40),rep("DGMM_comp3",40))
    
    plot(x=rep(data$mz,5),y=c(data$Average,data$Per_Pixel,data$DGMM_comp1,data$DGMM_comp2,data$DGMM_comp3),
         pch=16,cex=1, col=group,xlab="m/z", ylab="P value", main="P value" )
    abline(h=0.05, col="blue", lty=2)
    legend("topleft", legend=levels(group), pch=16, col=unique(group))
    
    
    
    
    p1<-ggplot(data, aes(mz, y = value)) + 
      geom_point(aes(y = Average, col = "Average"),size=2) + 
      geom_point(aes(y = Per_Pixel, col = "Per Pixel"),size=2)+
      geom_hline(yintercept=0.05, linetype="dashed", 
                 color = "blue", size=1)+
      ggtitle("P value") +
      xlab("m/z") + ylab("P value")
    
    return(list(av_p,pixel_p,dgmm_p))
  }
  
  if (compare=="ROC")
  {
    library("pROC")
    #########average
    
    prediction<-spectra(msset)
    category<-as.factor(msset$diagnosis)
    auc_av<-rep(0.5,length(mz(msset)))
    for(i in 1:40)
    {
      roc_obj <- roc(category, prediction[i,])
      auc_av[i]<-auc(roc_obj)
    }
    
    which(auc_av>sort(auc_av,decreasing=T)[10])
    
    ##############dgmm
    label<-list()
    for (f in 1:length(mz(msset)))
    {
      
      for ( i in 1:length(unique(msset$sample)))
      {
        dgmm<-fit_DGMM(msset=msset[,msset$sample=unique(msset$sample)[i]], f=f, k=0,out=3, kmax=3)
        label[[f]][i,]<-apply(dgmm$y,1, function (x) which(x==max(x)))
      }
    }
    ############match each mixture if mu is not ranked increasingly##########
    
    
    
    ########################################################################
    auc_dgmm<-matrix(0.5,ncol=kmax,nrow=length(mz(msset)))
    for (f in 1:length(mz(msset)))
    {
      k=length(unique(label[[f]]))
      print(k)
      for (i in 1:k)
      {
        prediction<-spectra(simImage)[f,label[[f]]==i]
        category<-as.factor(simImage$diagnosis[label[[f]]==i])
      }
      roc_obj <- roc(category, prediction)
      auc_dgmm[f,i]<-auc(roc_obj)
    }
    which(auc_dgmm>sort(auc_dgmm,decreasing=T)[10])
  }
   return(list(auc_av,auc_dgmm))
}