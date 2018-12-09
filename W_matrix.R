############single sample
W_matrix<-function(msset,sp_ratio=2,radius=1,sigma=1e4)
{
  ##################weight matrix
  radius<-1 
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
        w[i,j]<-exp(-((coord(rccnorep)$x[i]-coord(rccnorep)$x[j])^2+sp_ratio*(coord(rccnorep)$y[i]-coord(rccnorep)$y[j])^2))*exp(-t(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])%*%(spectra(rccnorep)[,j]-spectra(rccnorep)[,i])/sigma)
      }
    }
  }
  return(w)
}

