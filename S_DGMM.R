
####################step3
DGMM<-function(msset=msset,k=k,f=f,sp_ratio=4,step=1e5,)
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
        w[i,j]<-exp(-((coord(rccnorep)$x[i]-coord(rccnorep)$x[j])^2+sp_ratio*(coord(rccnorep)$y[i]-coord(rccnorep)$y[j])^2))
        *exp(-t(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])%*%(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])/10000)
      }
    }
  }
  
  ###################################fit DGMM
  Annl=1
  tt<-1
  int<-spectra(msset)[f,]
  gmm<-densityMclust(int,G=3,modelNames="V")
  x<-int
  N=length(x)
  k<-gmm$G
  g<-k
  #mu=c(99,130);
  #sigma=c(147,210);
  mu<-gmm$parameters$mean
  #mu<-c(30,60,90)
  sigma<-gmm$parameters$variance$sigmasq
  #sigma<-c(50,80,100)
  alpha=rep(1,g);
  beta=1;
  eta<-min(mu)/step
  dmu<-rep(1,g)
  dsg<-rep(1,g)
  dalpha<-rep(1,g)
  ##calculate priors
  y<-matrix(0, nrow=N, ncol=K)
  PI<-matrix(1/K, nrow=N, ncol=K)
  logPI<-matrix(1/K, nrow=N, ncol=K)
  px<-matrix(0, nrow=N, ncol=K)
  logpx<-matrix(0, nrow=N, ncol=K)
  iteration=100
  mutrace<-matrix(0,ncol=K,nrow=iteration)
  sigtrace<-matrix(0,ncol=K,nrow=iteration)
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
  y[y[,1]==0,1]<-rep(1e-200, length(y[y[,1]==0,1]))
  y[y[,2]==0,2]<-rep(1e-200, length(y[y[,2]==0,2]))
  
  for (i in 1:iteration)
  {
    
    
    
    ybar<-w%*%y/rowSums(w)
    # ybar[ybar[,1]==0,1]<-rep(1e-100, length(ybar[ybar[,1]==0,1]))
    # ybar[ybar[,2]==0,2]<-rep(1e-100, length(ybar[ybar[,2]==0,2]))
    ybar[ybar==0]<-1e-100
    loglik[i]<--sum(log(rowSums(t(t((ybar)^beta)*alpha^2)/rowSums(t(t((ybar)^beta)*alpha^2))*px)))
    logybar<-log(ybar)
    for ( j in 1:K)
    {
      logPI[,j]<-2*log(abs(alpha[j]))+beta*logybar[,j]
    }
    logPI<-logPI-rowMin(logPI)
    #  logPI<-pmin(logPI,100)
    #  logPI<-exp(logPI)
    #  PI<-logPI/rowSums(logPI)
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
    #  y[y[,1]==0,1]<-rep(1e-200, length(y[y[,1]==0,1]))
    #  y[y[,2]==0,2]<-rep(1e-200, length(y[y[,2]==0,2]))
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
  msset$dgmm<-apply(y,1, function (x) which(x==max(x)))
  image(msset, formula = dgmm~x*y,asp=5,colorkey=F)
  return(list(mu,sigma,alpha,beta,mutrace,sigtrace,alphatrace,betatrace,loglik,y))
}