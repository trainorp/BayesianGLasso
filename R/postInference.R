#' @export
postInference<-function(x) UseMethod("postInference")

#' @export
postInference.BayesianGLasso<-function(x)
{
  return(x$lambdas)
}

