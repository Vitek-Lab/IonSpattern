####################step 2




K_DGMM<-function(msset=msset,gmm=gmm,f=f,k=k,step=1e7,initialization="km",r_max=3)
{
  kr<-rep(0,r_max)
  for (radius in 1:r_max)
  {
    ##################neighboring matrix
    coords<-coord(msset)
    w <- apply(coords, 1, function(pt)
      (abs(as.numeric(pt["y"]) - as.numeric(coords[,'y'])) <= radius) & (abs(as.numeric(pt["x"]) - as.numeric(coords[,"x"])) <= radius))
    diag(w) <- F
    w<-apply(w,1,as.numeric)
    
    
    rmlist<-which(rowSums(w)==0)
    if (length(rmlist)!=0)
    {
      w<-w[-rmlist,-rmlist]
    }

    ###################################fit DGMM candidate
    ##########ion intensities
    int<-spectra(msset)[f,]
    if (length(rmlist)!=0)
    {
      int<-int[-rmlist]
    }
    x<-int
    #######number of pixels
    N=length(x)

    K<-k
    g<-k
    ###############initialize using k-means
    
    if (initialization=="km")
    {
      km<-kmeans(int,centers =k)
      mu<-km$centers
      sigma<-(mu*0.2)^2
      
    }
    
    ##############initialize using Gaussian Mixture Model
    
    if (initialization=="gmm")
    {
      if (gmm$G != k)
      {
        mu[1]<-0.1
        sigma[1]<-0.004
        mu[2:k]<-gmm$parameters$mean
        sigma[2:k]<-gmm$parameters$variance$sigmasq
      } else
      {
        mu<-gmm$parameters$mean
        sigma<-gmm$parameters$variance$sigmasq
      }
      
    }
    
    ############initialize alpha in Dirichlet process
    alpha=rep(1,g);
    ############initialize beta in PI
    beta=1;
    
    ###########step size
    eta<-min(mu)/step
    ##########differentials of mu, sigma and alpha
    dmu<-rep(1,g)
    dsg<-rep(1,g)
    dalpha<-rep(1,g)
    #########posterior probability
    y<-matrix(0, nrow=N, ncol=K)
    #########prior probability
    PI<-matrix(1/K, nrow=N, ncol=K)
    logPI<-matrix(1/K, nrow=N, ncol=K)
    #########P(x|mu, sigma)
    px<-matrix(0, nrow=N, ncol=K)
    logpx<-matrix(0, nrow=N, ncol=K)
    iteration=100
    #########trace
    mutrace<-matrix(0,ncol=K,nrow=iteration)
    sigtrace<-matrix(0,ncol=K,nrow=iteration)
    alphatrace<-matrix(0,ncol=K,nrow=iteration)
    betatrace<-matrix(0,ncol=1,nrow=iteration)
    #########negative loglikelihood
    loglik<-rep(0,iteration)
    #########initialize P(x|mu,sigma)
    for (j in 1:K)
    {
      px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
    }

    for (j in 1:K)
    {
      logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
    }
    ######### initialize posterior probability
    y<-px*PI/rowSums(px*PI)
    y[is.na(y)==TRUE]<-1/k
    y[y==0]<-1e-200
    
    for (i in 1:iteration)
    {
      
      
      ############average posterior probability
      ybar<-w%*%y/rowSums(w)
      
      ybar[ybar==0]<-1e-100
      
      #############negative loglikelihod
      loglik[i]<--sum(log(rowSums(t(t((ybar)^beta)*alpha^2)/rowSums(t(t((ybar)^beta)*alpha^2))*px)))
      logybar<-log(ybar)
      for ( j in 1:K)
      {
        logPI[,j]<-2*log(abs(alpha[j]))+beta*logybar[,j]
      }
      logPI<-logPI-rowMin(logPI)
      
      for ( j in 1:K)
      {
        PI[,j]<-alpha[j]^2*ybar[,j]^beta
      }
      PI[PI==Inf]<-1e100
      PI<-PI/rowSums(PI)
      
      PI[PI==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
      }
      ##logp(x|mu, sigma)
      for (j in 1:K)
      {
        logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
      }
      ##posterior
      
      
      y<-px*PI/rowSums(px*PI)
      y[is.na(y)==TRUE]<-1/k
      y[y==0]<-1e-100
      for ( j in 1:K)
      {
        dmu[j]<-sum(y[,j]*1/sigma[j]*(mu[j]-x))
        dsg[j]<-1/2*sum(y[,j]*(1/sigma[j]-(x-mu[j])^2/(sigma[j]^2)))
        dalpha_p<-y*(ybar[,j])^beta/rowSums(t(t((ybar)^beta)*alpha^2))
        dalpha_p[is.na(dalpha_p)]<-1
        dalpha[j]<--sum(2*y[,j]/alpha[j])+2*alpha[j]*sum(dalpha_p)
        
      }
      
      dbeta=sum(y*(-log(ybar)+rowSums(t(t((ybar)^beta)*alpha^2)*log(ybar))/rowSums(t(t((ybar)^beta)*alpha^2))))
      mu<-mu-eta*dmu
      sigma<-sigma-eta*dsg
      sigma[sigma<=0]<-0.006327605
      alpha<-alpha-eta*dalpha
      beta<-beta-eta*dbeta
      beta<-max(beta,0)
      beta<-min(beta,10)
      mutrace[i,]<-mu
      sigtrace[i,]<-sigma
      alphatrace[i,]<-alpha
      betatrace[i,]<-beta
      
    }
    xx<-rep(1,ncol(msset))
    if (length(rmlist)!=0)
    {
      xx[-rmlist]<-apply(ybar,1, function (x) which(x==max(x))) 
    } else
    {
      xx<-apply(ybar,1, function (x) which(x==max(x)))
    }

    msset$dgmm<-xx
#    image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F)
    kprim<-length(unique(msset$dgmm))
    L<-unique(msset$dgmm)
    for (i in L)
    {
      if (length(msset$dgmm[msset$dgmm==i])/ncol(msset)<0.01)
      {
        kprim<-kprim-1
        msset$dgmm[msset$dgmm==i]<-NA
      }
      
    }
    print(kprim)
    print(min(loglik))
    print(mu)
    seg<-unique(msset$dgmm[!is.na(msset$dgmm)])
    seg2<-seg
    for (i in 1:length(seg[!is.na(seg)]))
    {
      j<-i+1;
      while(j<length(seg[!is.na(seg)]))
      {
        if ((mu[seg[i]]-mu[seg[j]])/sqrt(sigma[seg[i]]/2+sigma[seg[j]]/2)<1)
        {
          seg2[j]<-NA
        }
        j=j+1
      }
    }
    kr[radius]<-length(seg2[!is.na(seg)])
  }
  kr<-Mode(kr)
  print(kr)
  return(kr)
}


