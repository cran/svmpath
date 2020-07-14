"InitsvmPath" <-
  function(Rmat, cvec, const)
{
###uses package kernlab
    n <- length(cvec)
    sol <- ipop(c=cvec,H=Rmat,A=matrix(1,1,n),b=const,l=rep(0,n),u=rep(1,n),r=0,margin=1e-5)
    list(alpha=sol@primal)
  }

