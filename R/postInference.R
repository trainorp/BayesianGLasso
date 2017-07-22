#' Posterior inference 
#'
#' Posterior distributions
#' 
#' @export
postInference<-function(x) UseMethod("postInference")

#' @export
postInference.BayesianGLasso<-function(x)
{
  return(x$lambdas)
}

