####################Matching
####simple
###########reference is ion image
label_ref<-apply(y,1, function (x) which(x==max(x)))
label_comp<-apply(y,1, function (x) which(x==max(x)))
MATCH<-function(msset=msset, label_ref=label_ref, label_comp=label_comp)
{
  n1<-length(unique(label_comp))
  n2<-length(unique(label_ref))
  match_score<-matrix(0,ncol=n2,nrow=n1)
  for ( i in unique(label_comp))
  {
    l1<-which(label_comp==i)
    for ( j in unique(label_ref))
    {
      l2<-which(label_ref==j)
      match_score[i,j]<-length(intersect(l1,l2))/(length(l2)/2+length(l1)/2)


      if (abs(match_score[i,j]-1)<0.05)
      {
        print(paste0("match",i,j))
      }
    }
  }
  print(match_score)
}


msset$dgmm<-(label_ref!=2)&(label_comp==1)
image(msset, formula = dgmm~x*y,asp=3,colorkey=F)

msset$dgmm<-label_ref
image(msset, formula = dgmm~x*y,asp=2,colorkey=F)