#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// using namespace arma; 

double rInvG(double mu, double lambda){
  double val = 0;
  double z,y,x,u;
  z = rnorm(1)[0];
  y = z*z;
  x = mu + 0.5*mu*mu*y/lambda - 0.5*(mu/lambda) * sqrt(4*mu*lambda*y + mu*mu*y*y);
  u = runif(1)[0];
  if(u <= mu/(mu+x)){
    val = x;
  }else{
    val = mu*mu/x;
  }
  return(val);
}

arma::vec upperTri(arma::mat a){
  arma::mat b = trimatu(a);
  arma::mat b2 = trimatl(a);
  arma::uvec c = find(b);
  arma::uvec c2 = find(b2);
  NumericVector d = wrap(c);
  NumericVector d2 = wrap(c2);
  NumericVector e = setdiff(d,d2);
  e=e.sort();
  arma::uvec f = as<arma::uvec>(e);
  arma::vec g = a.elem(f);
  return g; 
}

IntegerMatrix permFun(int p){
  IntegerVector permInt = seq_len(p) - 1;
  IntegerMatrix perms(p-1,p);
  IntegerVector permInt2(p);
  for(int i = 0; i < p; i++)
  {
    permInt2 = permInt[permInt != i];
    perms(_,i) = permInt2;
  }
  return perms;
}

// [[Rcpp::export]]
List bgl(int n, int iters, double lambdaPriorb, double lambdaPosta,
         arma::mat S, arma::mat Sigma, arma::mat Omega, int keepLambdas){
  // ProfilerStart("myprof.log");

  int p = Omega.n_cols;
  IntegerMatrix idk = permFun(p);
  List OmegaList(iters);
  
  NumericVector lambdaList;
  if(keepLambdas==1){
    lambdaList = NumericVector(iters);
  }
  else{
    lambdaList = NumericVector(1);
  }

  for(int iter=0; iter < iters; ++iter){

    // Sample lambda:
    arma::vec Omega2 = vectorise(Omega);
    double lambdaPostb = lambdaPriorb + sum(abs(Omega2)) / 2.0;
    double lambda = R::rgamma(lambdaPosta,1.0 / lambdaPostb);
    
    if(keepLambdas==1){
      lambdaList(iter) = lambda;
    }

    // Mu prime:
    arma::vec OmegaTemp = abs(upperTri(Omega));
    // arma::vec OmegaTemp = abs(Omega.elem(find(trimatu(Omega,1))));
    arma::vec mup = lambda/OmegaTemp;
    /* double mupThresh = pow(10.0,12.0);
    arma::uvec mupReplace = find(mup > mupThresh);
    mup.elem(mupReplace).fill(mupThresh); */

    // Sample tau:
    arma::vec rIG(p*(p-1)/2);
    for(unsigned int i = 0; i < mup.n_elem; i++){
      rIG(i) = rInvG(mup[i],lambda*lambda);
    }
    arma::vec tauVec = 1.0 / rIG;
    arma::mat tau(Omega.n_rows,Omega.n_cols,arma::fill::zeros);
    int k = 0;
    for(unsigned int j = 0; j < tau.n_cols; j++){
      for(unsigned int i = 0; i < j; i++){
        tau(i,j) = tauVec(k);
        tau(j,i) = tauVec(k);
        k++;
      }
    }

    IntegerVector tauIPerms(p);
    arma::uvec tauIPerms2(p), tauII(1);
    arma::vec Sigma12(p-1), rnorm1(p-1), beta(p-1);
    arma::mat tauI(p-1,1), Sigma11(p-1,p-1), Omega11inv(p-1,p-1), Ci(p-1,p-1), CiChol(p-1,p-1), mui(p-1,1);
    arma::mat OmegaInvTemp(p-1,1);
    double gamm;
    arma::mat gamm2;
    for(int i = 0; i < p; i++){
      tauIPerms = idk(_,i);
      tauIPerms2 = as<arma::uvec>(tauIPerms);
      tauII(0) = i;
      tauI = tau.submat(tauIPerms2, tauII);
      Sigma11 = Sigma.submat(tauIPerms2, tauIPerms2);
      Sigma12 = Sigma.submat(tauIPerms2, tauII);
      Omega11inv = Sigma11 - ((Sigma12 * Sigma12.t()) / Sigma(i,i));
      Ci = (S(i,i) + lambda) * Omega11inv + diagmat(1 / tauI);
      CiChol = arma::chol(Ci);
      mui = arma::solve(-Ci, S.submat(tauIPerms2, tauII));
      rnorm1 = rnorm(p-1);
      beta = mui + arma::solve(CiChol, rnorm1);
      Omega.submat(tauIPerms2, tauII) = beta;
      Omega.submat(tauII, tauIPerms2) = beta.t();
      gamm = rgamma(1,n/2+1,2/(S(i,i) + lambda))[0];
      gamm2 = gamm + beta.t() * Omega11inv * beta;
      Omega(i,i) = gamm2(0,0);
      OmegaInvTemp = Omega11inv * beta;
      Sigma.submat(tauIPerms2, tauIPerms2) = Omega11inv + (OmegaInvTemp * OmegaInvTemp.t())/gamm;
      Sigma.submat(tauIPerms2, tauII) = -OmegaInvTemp/gamm;
      Sigma.submat(tauII,tauIPerms2) = -OmegaInvTemp.t()/gamm;
      Sigma(i,i) = 1/gamm;
    }
    OmegaList(iter) = Omega;
    Rcpp::Rcout << "iter = " << iter + 1 << std::endl;

  }
  // ProfilerStop();

  return List::create(Named("lambdas") = lambdaList,
                      Named("Omegas") = OmegaList);
}

// [[Rcpp::export]]
List bAdgl(int n, int iters, double gammaPriors, double gammaPriort,
              double lambdaii, arma::mat S, arma::mat Sigma,
              arma::mat Omega, bool priorLogical, arma::mat priorHyper,
              int keepLambdas){
  // ProfilerStart("myprof.log");

  int p = Omega.n_cols;
  IntegerMatrix idk = permFun(p);

  List OmegaList(iters);

  List LambdaList;
  if(keepLambdas==1){
    LambdaList = List(iters);
  }
  else{
    LambdaList = List(1);
  }

  for(int iter=0; iter < iters; ++iter){
    // Sample lambdas:
    arma::vec Omega2 = vectorise(Omega);
    arma::vec OmegaTemp = abs(upperTri(Omega));
    // arma::vec OmegaTemp = abs(Omega.elem(find(trimatu(Omega,1))));

    // Gamma distirbution conditional parameters:
    arma::vec tt = OmegaTemp + gammaPriort;

    // If informative prior:
    if(priorLogical){
      arma::vec priorHyperVec = vectorise(priorHyper);
      priorHyperVec = upperTri(priorHyper);
      // priorHyperVec = priorHyperVec.elem(find(trimatu(priorHyper,1)));
      tt = tt + priorHyperVec;
    }

    arma::vec s(tt.n_elem);
    s.fill(gammaPriors + 1);

    // Sample lambda:
    arma::vec lambda(tt.n_elem);
    for(unsigned int ii = 0; ii<tt.n_elem; ++ii){
      lambda(ii) = R::rgamma(s(ii), 1.0/tt(ii)); // Change to Rcpp sugar version
    }

    // Mu prime:
    arma::vec mup = lambda/OmegaTemp;
    double mupThresh = pow(10.0,12.0);
    arma::uvec mupReplace = find(mup > mupThresh);
    mup.elem(mupReplace).fill(mupThresh);

    // Sample tau:
    arma::vec rIG(p*(p-1)/2);
    for(unsigned int i = 0; i < mup.n_elem; i++){
      rIG(i) = rInvG(mup(i),pow(lambda(i),2));
    }
    arma::vec tauVec = 1.0 / rIG;
    arma::mat tau(Omega.n_rows,Omega.n_cols,arma::fill::zeros);
    int k = 0;
    for(unsigned int j = 0; j < tau.n_cols; j++){
      for(unsigned int i = 0; i < j; i++){
        tau(i,j) = tauVec(k);
        tau(j,i) = tauVec(k);
        k++;
      }
    }

    IntegerVector tauIPerms(p);
    arma::uvec tauIPerms2(p), tauII(1);
    arma::vec Sigma12(p-1), rnorm1(p-1), beta(p-1);
    arma::mat tauI(p-1,1), Sigma11(p-1,p-1), Omega11inv(p-1,p-1), Ci(p-1,p-1), CiChol(p-1,p-1), mui(p-1,1);
    arma::mat OmegaInvTemp(p-1,1);
    double gamm;
    arma::mat gamm2;
    for(int i = 0; i < p; i++){
      tauIPerms = idk(_,i);
      tauIPerms2 = as<arma::uvec>(tauIPerms);
      tauII(0) = i;
      tauI = tau.submat(tauIPerms2, tauII);
      Sigma11 = Sigma.submat(tauIPerms2, tauIPerms2);
      Sigma12 = Sigma.submat(tauIPerms2, tauII);
      Omega11inv = Sigma11 - ((Sigma12 * Sigma12.t()) / Sigma(i,i));
      Ci = (S(i,i) + lambdaii) * Omega11inv + diagmat(1 / tauI);
      
      CiChol = arma::chol(Ci);
      mui = arma::solve(-Ci, S.submat(tauIPerms2, tauII));
      rnorm1 = rnorm(p-1);
      beta = mui + arma::solve(CiChol, rnorm1);
      Omega.submat(tauIPerms2, tauII) = beta;
      Omega.submat(tauII, tauIPerms2) = beta.t();
      gamm = rgamma(1,n/2+1,2/(S(i,i) + lambdaii))[0];
      gamm2 = gamm + beta.t() * Omega11inv * beta;
      Omega(i,i) = gamm2(0,0);
      OmegaInvTemp = Omega11inv * beta;
      Sigma.submat(tauIPerms2, tauIPerms2) = Omega11inv + (OmegaInvTemp * OmegaInvTemp.t())/gamm;
      Sigma.submat(tauIPerms2, tauII) = -OmegaInvTemp/gamm;
      Sigma.submat(tauII,tauIPerms2) = -OmegaInvTemp.t()/gamm;
      Sigma(i,i) = 1/gamm;
    }
    OmegaList(iter) = Omega;
    
    arma::mat lambdaMat(Omega.n_rows,Omega.n_cols,arma::fill::zeros);
    int kk = 0;
    for(unsigned int j = 0; j < lambdaMat.n_cols; j++){
      for(unsigned int i = 0; i < j; i++){
        lambdaMat(i,j) = lambda(kk);
        kk++;
      }
    }
    if(keepLambdas==1){
      LambdaList(iter) = lambdaMat;
    }
    
    double detOmega = arma::det(Omega);
    Rcpp::Rcout << "iter = " << iter + 1 << " det(Omega) = " << detOmega << std::endl;
  }

  return List::create(Named("lambdas") = LambdaList,
                      Named("Omegas") = OmegaList);
}

