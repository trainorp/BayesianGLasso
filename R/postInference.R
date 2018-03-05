#' Posterior inference 
#'
#' Posterior distributions of [LOH]
#' 
#' @param x Bayesian Graphical Lasso object as returned by the function blockGLasso
#' @export
posteriorInference<-function(x,...) UseMethod("posteriorInference")

#' @export
posteriorInference.BayesianGLasso<-function(x,parameters="Omega",alpha=.05)
{
  # Check parameters
  if(! parameters %in% c("Omega","Sigma"))
  {
    stop("Parameters argument must be Omega (concentration matrix) or 
         Sigma (covariance matrix)")
  }

  stats<-c("mean","sd","length","median")

  # Discard burn in:
  x2<-x
  x2$Sigmas<-x2$Sigmas[(x2$burnIn+1):length(x2$Sigmas)]
  x2$Omegas<-x2$Omegas[(x2$burnIn+1):length(x2$Omegas)]
  if(!is.null(x2$lambdas)) x2$lambdas<-x2$lambdas[(x2$burnIn+1):length(x2$lambdas)]
  
  for(stat in stats)
  {
    statMatrix<-matrix(NA,nrow=nrow(x2$Omegas[[1]]),ncol=ncol(x2$Omegas[[1]]))
    for(i in 1:nrow(statMatrix))
    {
      for(j in 1:ncol(statMatrix))
      {
        statMatrix[i,j]<-eval(call(stat,sapply(x2$Omegas,FUN=function(y) y[i,j])))
      }
    }
    assign(paste0("stat",stat),statMatrix)
  }
  
  lq<-function(z) quantile(z,probs=alpha/2)
  uq<-function(z) quantile(z,probs=1-(alpha/2))
  lqMatrix<-uqMatrix<-matrix(NA,nrow=nrow(x2$Omegas[[1]]),ncol=ncol(x2$Omegas[[1]]))
  for(i in 1:nrow(lqMatrix))
  {
    for(j in 1:ncol(lqMatrix))
    {
      lqMatrix[i,j]<-eval(call("lq",sapply(x2$Omegas,FUN=function(y) y[i,j])))
      uqMatrix[i,j]<-eval(call("uq",sapply(x2$Omegas,FUN=function(y) y[i,j])))
    }
  }
  
  list(posteriorMean=statmean,posteriorSd=statsd,posteriorMedian=statmedian,
       lowerCI=lqMatrix,upperCI=uqMatrix)
}

