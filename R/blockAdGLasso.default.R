#' @export
blockAdGLasso.default<-function(X,iterations=2000,burnIn=1000,adaptiveType=c("norm","priorHyper"),
                                priorHyper=NULL,gammaPriors=1,gammaPriort=1,
                                lambdaii=1,illStart=c("identity","glasso"),rho=.1,
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

  # Adaptive type:
  adaptiveType<-match.arg(adaptiveType)
  if(adaptiveType=="priorHyper"){
    if(is.null(priorHyper)) stop("Must specify matrix of prior hyperparameters")
    if(class(priorHyper)!="matrix") stop("Prior hyperparameters must be provided as a matrix")
    if(!all(dim(Sigma)==dim(priorHyper))) stop("Dimesion of hyperparameters does not equal dimension
                                           of the concentration matrix")
    priorLogical<-TRUE
  }else{
    priorLogical<-FALSE
    priorHyper<-matrix(1,nrow=p,ncol=p)
  }

  # Concentration matrix and it's dimension:
  if(rcond(Sigma)<.Machine$double.eps){
    if(illStart=="identity"){
      Omega<-diag(nrow(Sigma))+1/(p**2)
    }
    else{
      Omega<-glasso::glasso(cov(X),rho=rho)$wi+1/(p**2)
    }
  }else{
    Omega<-MASS::ginv(Sigma)
  }
  rownames(Omega)<-rownames(Sigma)
  colnames(Omega)<-colnames(Sigma)

  bglObj<-bAdgl(n=n, iters=totIter, gammaPriors=gammaPriors,
                gammaPriort=gammaPriort, lambdaii=lambdaii, S=S, Sigma=Sigma,
                Omega=Omega, priorLogical=priorLogical, priorHyper=priorHyper)

  bglObj<-c(bglObj,burnIn=burnIn)
  class(bglObj)<-"BayesianGLasso"
  return(bglObj)
}
