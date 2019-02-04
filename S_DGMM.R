
####################step3
DGMM<-function(msset=msset,w=w,k=k,f=f,sp_ratio=4,step=1e5, iteration=1000,Annl=0, sig=0.2,initialization="Km")
{

  ############# remove isolated pixels
  
  rmlist<-which(rowSums(w)==0)
  
  
  if (length(rmlist)!=0)
  {
    w<-w[-rmlist,-rmlist]
  }
  
  ###################################fit DGMM
  int<-spectra(msset)[f,]
  if (length(rmlist)!=0)
  {
    int<-int[-rmlist]
  }
  #Annl=0
  tt<-1
  
  x<-int
  N=length(x)
  ###############initialize using k-means
  
  if (initialization=="km")
  {
    km<-kmeans(int,centers =k)
    mu<-sort(km$centers, decreasing = F)
    sigma<-(mu*sig)^2
    sigma[sigma<0.0006]<-0.0006
    
  }
  
  ##############initialize using Gaussian Mixture Model
  
  if (initialization=="gmm")
  {
    gmm<-densityMclust(int,G=k,modelNames="V")
    mu<-gmm$parameters$mean
    sigma<-gmm$parameters$variance$sigmasq
  }
  
  ############initialize alpha in Dirichlet process
  
  alpha<-rep(1,k)
  ############initialize beta in PI
  beta=1;
  K<-k
  #########step size
  eta<-min(mu)/step
  ###########differentials of mu, sigma and alpha
  dmu<-rep(1,k)
  dsg<-rep(1,k)
  dalpha<-rep(1,k)
  ###########posterior probability
  y<-matrix(0, nrow=N, ncol=K)
  #########prior probability
  PI<-matrix(1/K, nrow=N, ncol=K)
  logPI<-matrix(1/K, nrow=N, ncol=K)
  #########P(x|mu, sigma)
  px<-matrix(0, nrow=N, ncol=K)
  logpx<-matrix(0, nrow=N, ncol=K)
  #iteration=100
  #########trace
  mutrace<-matrix(0,ncol=k,nrow=iteration)
  sigtrace<-matrix(0,ncol=k,nrow=iteration)
  alphatrace<-matrix(0,ncol=K,nrow=iteration)
  betatrace<-matrix(0,ncol=1,nrow=iteration)
  #########negative loglikelihood
  loglik<-rep(0,iteration)
  PI_p<-matrix(0, nrow=N, ncol=K)
  px_p<-matrix(0, nrow=N, ncol=K)
  #########initialize P(x|mu,sigma)
  for (j in 1:K)
  {
    px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
  }
  ##logp(x|\theta)
  for (j in 1:K)
  {
    logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
  }
  ######## initialize posterior probability
  
  y<-px*PI/rowSums(px*PI)
  y[is.na(y)==TRUE]<-1/k
  ##########handling data out of storage range
  y[y==0]<-1e-100
  
  
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
    ##logp(x|mu, sigma)1
    for (j in 1:K)
    {
      logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
    }
    ##posterior
    
    
    y<-px*PI/rowSums(px*PI)
    y[is.na(y)==TRUE]<-1/k
    y[y==0]<-1e-100
    
    ###################calculate differential
    for ( j in 1:K)
    {
      dmu[j]<-sum(y[,j]*1/sigma[j]*(mu[j]-x))
      dsg[j]<-1/2*sum(y[,j]*(1/sigma[j]-(x-mu[j])^2/(sigma[j]^2)))
      dalpha_p<-y*(ybar[,j])^beta/rowSums(t(t((ybar)^beta)*alpha^2))
      dalpha_p[is.na(dalpha_p)]<-1
      dalpha[j]<--sum(2*y[,j]/alpha[j])+2*alpha[j]*sum(dalpha_p)
      
    }
    
    dbeta=sum(y*(-log(ybar)+rowSums(t(t((ybar)^beta)*alpha^2)*log(ybar))/rowSums(t(t((ybar)^beta)*alpha^2))))
    ##################updata parameters
    mu<-mu-eta*dmu
    for (mu_ind in 1:k)
    {
      mu[mu_ind]<-max(mu[mu_ind],0)
    }
    sigma<-sigma-eta*dsg
    sigma[sigma<=0]<-0.006327605
    alpha<-alpha-eta*dalpha
    beta<-beta-eta*dbeta
    beta<-max(beta,0)
    beta<-min(beta,10)

    ############################simulated annealing
    if (Annl==TRUE)
    {
      ####randomly change first element of mu
      mu_p<-c(runif(1, min=max(mu[1]-0.2*mu[1]*tt,0), max=mu[1]+0.2*mu[1]*tt), mu[2:k])
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2
      }
      PI_p<-PI/rowSums(PI)

      PI_p[PI_p==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu_p[j])^2/2/sigma[j])
      }
      y_p<-PI_p*px_p/rowSums(PI_p*px_p)
      y_p[y_p==0]<-1e-100
      ybar_p<-w%*%y_p/rowSums(w)
      ybar_p[ybar_p==0]<-1e-100
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2*ybar_p[,j]^beta
      }
      PI_p<-PI/rowSums(PI)

      PI_p[PI_p==0]<-1e-100 
      ##p(x|mu, sigma)
      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu_p[j])^2/2/sigma[j])
      }
      
      loglik1<--sum(log(rowSums(t(t((ybar_p)^beta)*alpha^2)/rowSums(t(t((ybar_p)^beta)*alpha^2))*px_p)))
      #######################
      ybar_p<-w%*%y/rowSums(w)
      ybar_p[ybar_p==0]<-1e-100
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2*ybar[,j]^beta
      }
      PI_p<-PI/rowSums(PI)

      PI_p[PI_p==0]<-1e-100 

      for (j in 1:K)
      {
        px_p[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
      }
      
      loglik2<--sum(log(rowSums(t(t((ybar_p)^beta)*alpha^2)/rowSums(t(t((ybar_p)^beta)*alpha^2))*px_p)))
      
      
      if (loglik1<loglik2)
      {
        mu<-mu_p
        y<-y_p
      }
    }
    ###################################################################
    
    mutrace[i,]<-mu
    sigtrace[i,]<-sigma
    alphatrace[i,]<-alpha
    betatrace[i,]<-beta
    #########cooling
    tt<-tt-1/iteration


  }

  
  
  xx<-rep(1,ncol(msset))
  if (length(rmlist)!=0)
  {
    xx[-rmlist]<-apply(y,1, function (x) which(x==max(x))) 
  } else
  {
    xx<-apply(y,1, function (x) which(x==max(x)))
  }
  
  msset$dgmm<-xx
#  image(msset, formula = dgmm~x*y,asp=sp_ratio, main=paste0(f,"y"))
  
  xx<-rep(1,ncol(msset))
  if (length(rmlist)!=0)
  {
    xx[-rmlist]<-apply(ybar,1, function (x) which(x==max(x))) 
  } else
  {
    xx<-apply(y,1, function (x) which(x==max(x)))
  }
  
  msset$dgmm<-xx
#  image(msset, formula = dgmm~x*y,asp=sp_ratio, main=paste0(f,"ybar"))

  return(list(mu,sigma,alpha,beta,mutrace,sigtrace,alphatrace,betatrace,loglik,y))
}