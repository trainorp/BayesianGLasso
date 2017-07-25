#' Posterior inference 
#'
#' Posterior distributions
#' 
#' @export
posteriorCI<-function(x,...) UseMethod("posteriorCI")

#' @export
posteriorCI.BayesianGLasso<-function(x,parameters="Omega")
{
  # Check parameters
  if(! parameters %in% c("Omega","Sigma"))
  {
    stop("Parameters argument must be Omega (concentration matrix) or 
         Sigma (covariance matrix)")
  }
  
  # Discard burn in:
  x2<-x
  x2$Sigmas<-x2$Sigmas[(x2$burnIn+1):length(x2$Sigmas)]
  x2$Omegas<-x2$Omegas[(x2$burnIn+1):length(x2$Omegas)]
  if(!is.null(x2$lambdas)) x2$lambdas<-x2$lambdas[(x2$burnIn+1):length(x2$lambdas)]
  
  print(x$burnIn)
  return(x$lambdas)
}

