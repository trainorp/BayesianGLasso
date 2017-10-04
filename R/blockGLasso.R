#' Block Gibbs sampler for Bayesian Graphical Lasso
#'
#' Blockwise sampling from the conditional distribution of a permuted column/row
#' for simulating the posterior distribution for the concentration matrix specifying
#' a Gaussian Graphical Model
#' @param X Data matrix
#' @param iterations Length of Markov chain after burn-in
#' @param burnIn Number of burn-in iterations
#' @param adaptive Logical; Adaptive graphical lasso (TRUE) or regular (FALSE). Default is FALSE.
#' @param lambdaPriora Shrinkage parameter (lambda) gamma distribution shape hyperparameter
#' (Ignored if adaptive=TRUE)
#' @param lambdaPriorb Shrinkage parameter (lambda) gamma distribution scale hyperparameter
#' (Ignored if adaptive=TRUE)
#' @param adaptiveType Choose of adaptive type. Options are "norm" for norm of concentration
#' matrix based adaptivity and "priorHyper" for informative adaptivity
#' @param priorHyper Matrix of gamma scale hyper parameters (Ignord if adaptiveType="norm")
#' @param gammaPriors labmda_ij gamma distribution shape prior (Ignored if adaptive=FALSE)
#' @param gammaPriort lambda_ij gamma distribution rate prior (Ignored if adaptive=FALSE)
#' @param lambdaii lambda_ii hyperparameter (Ignored if adaptive=FALSE)
#' @param illStart Method for generating a positive definite estimate of the sample covariance matrix if sample covariance matrix is not semi-positive definite
#' @param rho Regularization parameter for the graphical lasso estimate of the sample covariance matrix (if illStart="glasso")
#' @param verbose logical; if TRUE return MCMC progress
#' @details Implements the block Gibbs sampler for the Bayesian graphical lasso
#' introduced in Wang (2012). Samples from the conditional distribution of a 
#' permuted column/row for simulating the posterior distribution for the concentration 
#' matrix specifying a Gaussian Graphical Model
#' @return 
#' \item{Sigma}{List of covariance matrices from the Markov chain}
#' \item{Omega}{List of concentration matrices from the Markov chains}
#' \item{Lambda}{Vector of simulated lambda parameters}
#' @author Patrick Trainor (University of Louisville)
#' @author Hao Wang
#' @references Wang, H. (2012). Bayesian graphical lasso models and efficient 
#' posterior computation. \emph{Bayesian Analysis, 7}(4). <doi:10.1214/12-BA729> .
#' @examples
#' \donttest{
#' # Generate true covariance matrix:
#' s<-.9**toeplitz(0:9)
#' # Generate multivariate normal distribution:
#' set.seed(5)
#' x<-MASS::mvrnorm(n=100,mu=rep(0,10),Sigma=s)
#' blockGLasso(X=x)
#' }
#' # Same example with short MCMC chain:
#' s<-.9**toeplitz(0:9)
#' set.seed(6)
#' x<-MASS::mvrnorm(n=100,mu=rep(0,10),Sigma=s)
#' blockGLasso(X=x,iterations=100,burnIn=100)
#' @export
blockGLasso<-function(X,iterations=2000,burnIn=1000,adaptive=FALSE,
          lambdaPriora=1,lambdaPriorb=1/10,
          adaptiveType=c("norm","priorHyper"), priorHyper=NULL,
          gammaPriors=1,gammaPriort=1,
          lambdaii=1,illStart=c("identity","glasso"),rho=.1,
          verbose=TRUE,...)
{
  if(adaptive)
  {
    UseMethod("blockAdGLasso")
  }else{
    UseMethod("blockGLasso")
  }
}
