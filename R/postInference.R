#' Posterior inference 
#'
#' Inferential statistics for the simulated posterior distributions of the covariance 
#' and concentration parameter matrices for Gaussian Graphical Models estimated
#' via the adaptive or non-adaptive Bayesian Graphical Lasso
#' 
#' @param x Bayesian Graphical Lasso object as returned by the function blockGLasso
#' @param parameter Options are "Omega" for the posterior distribution of the 
#' concentration matrix or "Sigma" for the covariance matrix
#' @return
#' \item{posteriorMean}{Matrix of mean values of the simulated posterior distribution}
#' \item{posteriorSd}{Matrix of standard deviation values of the simulated posterior distribution}
#' \item{posteriorMedian}{Matrix of median values of the simulated posterior distribution}
#' \item{lowerCI}{Lower limit of the Bayesian credible interval for each matrix parameter}
#' \item{upperCI}{Upper limit of the Bayesian credible interval for each matrix parameter}
#' @export
posteriorInference<-function(bglObj,parameter="Omega",
                             alpha=.05,...){
  UseMethod("posteriorInference")
}

#' @export
posteriorInference.BayesianGLasso<-function(bglObj,parameter="Omega",
                                            alpha=.05){
  # Check parameters
  if(! parameter %in% c("Omega","Sigma")){
    stop("Parameter argument must be Omega (concentration matrix) or 
         Sigma (covariance matrix)")
  }

  stats<-c("mean","sd","length","median")

  # Discard burn in:
  x2<-bglObj
  x2$Sigmas<-x2$Sigmas[(x2$burnIn+1):length(x2$Sigmas)]
  x2$Omegas<-x2$Omegas[(x2$burnIn+1):length(x2$Omegas)]
  if(!is.null(x2$lambdas)) x2$lambdas<-x2$lambdas[(x2$burnIn+1):length(x2$lambdas)]
  
  for(stat in stats){
    statMatrix<-matrix(NA,nrow=nrow(x2$Omegas[[1]]),ncol=ncol(x2$Omegas[[1]]))
    for(i in 1:nrow(statMatrix)){
      for(j in 1:ncol(statMatrix)){
        statMatrix[i,j]<-eval(call(stat,sapply(x2$Omegas,FUN=function(y) y[i,j])))
      }
    }
    assign(paste0("stat",stat),statMatrix)
  }
  
  lq<-function(z) quantile(z,probs=alpha/2)
  uq<-function(z) quantile(z,probs=1-(alpha/2))
  lqMatrix<-uqMatrix<-matrix(NA,nrow=nrow(x2$Omegas[[1]]),ncol=ncol(x2$Omegas[[1]]))
  for(i in 1:nrow(lqMatrix)){
    for(j in 1:ncol(lqMatrix)){
      lqMatrix[i,j]<-eval(call("lq",sapply(x2$Omegas,FUN=function(y) y[i,j])))
      uqMatrix[i,j]<-eval(call("uq",sapply(x2$Omegas,FUN=function(y) y[i,j])))
    }
  }
  
  list(parameter=parameter,alpha=alpha,posteriorMean=statmean,posteriorSd=statsd,
       posteriorMedian=statmedian,lowerCI=lqMatrix,upperCI=uqMatrix)
}

