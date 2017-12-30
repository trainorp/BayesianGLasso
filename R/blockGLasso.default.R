#' @export
blockGLasso.default<-function(X,iterations=2000,burnIn=1000,lambdaPriora=1,lambdaPriorb=1/10,
                              illStart=c("identity","glasso"),rho=.1,
                              verbose=TRUE,...){
  # Total iterations:
  totIter<-iterations+burnIn

  # Ill conditioned start:
  illStart<-match.arg(illStart)

  # Sum of product matrix, covariance matrix, n
  S<-t(X)%*%X
  n=nrow(X)
  Sigma=S/n
  p<-dim(Sigma)[1]

  # Concentration matrix and it's dimension:
  if(rcond(Sigma)<.Machine$double.eps){
    if(illStart=="identity"){
      Omega<-diag(nrow(Sigma))+1/(p**2)
    }
    else{
      Omega<-glasso::glasso(cov(X),rho=rho)$wi+1/(p**2)
    }
  }
  else{
    Omega<-MASS::ginv(Sigma)
  }
  rownames(Omega)<-rownames(Sigma)
  colnames(Omega)<-colnames(Sigma)

  # Gamma distirbution posterior parameter a:
  lambdaPosta<-(lambdaPriora+(p*(p+1)/2))

  bglObj<-bgl(n=n,iters=totIter,lambdaPriorb=lambdaPriorb,lambdaPosta=lambdaPosta,
              S=S,Sigma=Sigma,Omega=Omega)

  bglObj<-c(bglObj,burnIn=burnIn)
  class(bglObj)<-"BayesianGLasso"
  return(bglObj)
}
