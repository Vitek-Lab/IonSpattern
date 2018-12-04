
####################step3
DGMM<-function(msset=msset,k=k,f=f,sp_ratio=4,step=1e5,Annl=0, initialization="Km")
{
  ##################weight matrix
  coords<-coord(msset)
  w <- apply(coords, 1, function(pt)
    (abs(as.numeric(pt["y"]) - as.numeric(coords[,'y'])) <= radius) & (abs(as.numeric(pt["x"]) - as.numeric(coords[,"x"])) <= radius))
  diag(w) <- F
  w<-apply(w,1,as.numeric)
  
  #################similarity score
  rccnorep<-msset
  for (i in 1:nrow(w))
  {
    for (j in 1:nrow(w))
    {
      
      if (w[i,j]==1)
      {
        w[i,j]<-exp(-((coord(rccnorep)$x[i]-coord(rccnorep)$x[j])^2+sp_ratio*(coord(rccnorep)$y[i]-coord(rccnorep)$y[j])^2))*exp(-t(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])%*%(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])/10000)
      }
    }
  }
  
  
  rmlist<-which(rowSums(w)==0)
  w<-w[-rmlist,-rmlist]
  ###################################fit DGMM
  
  
  #Annl=0
  tt<-1
  int<-spectra(msset)[f,]
  
  int<-int[-rmlist]
  x<-int
  N=length(x)
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
    gmm<-densityMclust(int,G=k,modelNames="V")
    mu<-gmm$parameters$mean
    sigma<-gmm$parameters$variance$sigmasq
  }
  
 ############initialize alpha beta
 
  alpha<-rep(1,k)
  beta=1;
  
  #########step size
  eta<-min(mu)/step
  ###########differential
  dmu<-rep(1,k)
  dsg<-rep(1,k)
  dalpha<-rep(1,k)
  ##calculate priors
  y<-matrix(0, nrow=N, ncol=K)
  PI<-matrix(1/K, nrow=N, ncol=K)
  logPI<-matrix(1/K, nrow=N, ncol=K)
  px<-matrix(0, nrow=N, ncol=K)
  logpx<-matrix(0, nrow=N, ncol=K)
  iteration=100
  mutrace<-matrix(0,ncol=k,nrow=iteration)
  sigtrace<-matrix(0,ncol=k,nrow=iteration)
  alphatrace<-matrix(0,ncol=K,nrow=iteration)
  betatrace<-matrix(0,ncol=1,nrow=iteration)
  loglik<-rep(0,iteration)
  PI_p<-matrix(0, nrow=N, ncol=K)
  px_p<-matrix(0, nrow=N, ncol=K)
  ##p(x|\theta)
  for (j in 1:K)
  {
    px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
  }
  ##logp(x|\theta)
  for (j in 1:K)
  {
    logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
  }
  
  y<-px*PI/rowSums(px*PI)
  ##########handling data out of storage range
  y[y[,1]==0,1]<-rep(1e-200, length(y[y[,1]==0,1]))
  y[y[,2]==0,2]<-rep(1e-200, length(y[y[,2]==0,2]))
  
  for (i in 1:iteration)
  {
    ############ybar
    ybar<-w%*%y/rowSums(w)

    ybar[ybar==0]<-1e-100
    
    #############loglikelihod
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
    ##p(x|\theta)
    for (j in 1:K)
    {
      px[,j]<-1/(2* pi )^0.5*1/sigma[j]^0.5*exp(-(x-mu[j])^2/2/sigma[j])
    }
    ##logp(x|\theta)
    for (j in 1:K)
    {
      logpx[,j]<-log(1/(2* pi )^0.5*1/sigma[j]^0.5)-(x-mu[j])^2/2/sigma[j]
    }
    ##posterior
    
    
    y<-px*PI/rowSums(px*PI)

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
    sigma<-sigma-eta*dsg
    sigma[sigma<=0]<-0.006327605
    alpha<-alpha-eta*dalpha
    beta<-beta-eta*dbeta
    beta<-max(beta,0)
    beta<-min(beta,10)

    ############################simulated annealing
    if (Annl==TRUE)
    {
      #  beta_p<-runif(1, min=max(beta-0.5*beta*tt,0), max=beta+0.5*beta*tt)
      #  beta_p<-runif(1,min=0, max=10)
      mu_p<-c(runif(1, min=max(mu[1]-0.2*mu[1]*tt,0), max=mu[1]+0.2*mu[1]*tt), mu[2:k])
      for ( j in 1:K)
      {
        PI_p[,j]<-alpha[j]^2
      }
      PI_p<-PI/rowSums(PI)
      #  PI[PI[,1]==0,1]<-rep(1e-100, length(PI[PI[,1]==0,1]))
      #  PI[PI[,2]==0,2]<-rep(1e-100, length(PI[PI[,2]==0,2]))
      PI_p[PI_p==0]<-1e-100 
      ##p(x|\theta)
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
        PI_p[,j]<-alpha[j]^2*ybar[,j]^beta
      }
      PI_p<-PI/rowSums(PI)
      #  PI[PI[,1]==0,1]<-rep(1e-100, length(PI[PI[,1]==0,1]))
      #  PI[PI[,2]==0,2]<-rep(1e-100, length(PI[PI[,2]==0,2]))
      PI_p[PI_p==0]<-1e-100 
      ##p(x|\theta)
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
      #  PI[PI[,1]==0,1]<-rep(1e-100, length(PI[PI[,1]==0,1]))
      #  PI[PI[,2]==0,2]<-rep(1e-100, length(PI[PI[,2]==0,2]))
      PI_p[PI_p==0]<-1e-100 
      ##p(x|\theta)
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
  xx[-rmlist]<-apply(y,1, function (x) which(x==max(x)))
  msset$dgmm<-xx
  image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F)
  
  msset$dgmm<-y[,1]
  image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F,main=paste0(i))
  msset$dgmm<-apply(ybar,1, function (x) which(x==max(x)))
  image(msset, formula = dgmm~x*y,asp=sp_ratio,colorkey=F,main=paste0(i))

  return(list(mu,sigma,alpha,beta,mutrace,sigtrace,alphatrace,betatrace,loglik,y))
}