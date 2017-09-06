#' @export
blockAdGLasso.default<-function(X,iterations=2000,burnIn=1000,lambdaPriora=1,lambdaPriorb=1/10,
                                verbose=TRUE,...)
{
  # Total iterations:
  totIter<-iterations+burnIn 
  
  # Sum of product matrix, covariance matrix, n
  S<-t(X)%*%X
  n=nrow(X)
  Sigma=S/n
  
  # Concentration matrix and it's dimension:
  Omega<-MASS::ginv(Sigma)
  p<-dim(Omega)[1]
  
  # Indicator matrix and permutation matrix for looping through columns & rows ("blocks")
  indMat<-matrix(1:p**2,ncol=p,nrow=p)
  perms<-matrix(NA,nrow=p-1,ncol=p)
  permInt<-1:p
  for(i in 1:ncol(perms))
  {
    perms[,i]<-permInt[-i]
  }
  
  # Structures for storing each MCMC iteration:
  SigmaMatList<-OmegaMatList<-list()

  # Latent tau:
  tau<-matrix(NA,nrow=p,ncol=p)
  
  # lambda hyperparameters:
  s<-1e-2  
  tt<-1e-6
  lambdaii<-1
  
  # Main block sampling loop:
  tryCatch(
    {
      for(iter in 1:totIter)
      {
        print(iter)
        OmegaTemp<-Omega[lower.tri(Omega)]
        OmegaTemp<-abs(OmegaTemp)
        OmegaTemp<-ifelse(OmegaTemp<1e-6,1e-6,OmegaTemp)
        
        # Gamma distirbution posterior parameter s:
        s<-s+1
        
        # Gamma distirbution posterior parameter t:
        tt<-OmegaTemp+tt
        
        # Sample lambda (LOH):
        lambda<-sapply(tt,FUN=function(x) stats::rgamma(1,shape=s,scale=1/x))
        
        mup<-lambda/OmegaTemp
        mup<-ifelse(mup>10e12,10e12,mup)
        
        # Sample tau:
        rIG<-rep(NA,length(mup))
        for(ii in 1:length(mup))
        {
          rIG[ii]<-statmod::rinvgauss(n=1,mean=mup[ii],shape=lambda[ii]**2)
        }

        tau[upper.tri(tau)]<-1/rIG
        tau[lower.tri(tau)]<-t(tau)[lower.tri(t(tau))]
        
        # Sample from conditional distribution by column:
        for(i in 1:p)
        {
          cat("i is ", i)
          tauI<-tau[perms[,i],i]
          Sigma11<-Sigma[perms[,i],perms[,i]]
          Sigma12<-Sigma[perms[,i],i]
          
          Omega11inv<-Sigma11-Sigma12%*%t(Sigma12)/Sigma[i,i]
          Ci<-(S[i,i]+lambdaii)*Omega11inv+diag(1/tauI)
          
          CiChol<-chol(Ci)
          mui<-solve(-Ci,S[perms[,i],i])
          
          # Sampling:
          rnorm1<-stats::rnorm(p-1)
          beta<-mui+solve(CiChol,rnorm1)
          
          # Replacing omega entries
          Omega[perms[,i],i]<-beta
          Omega[i,perms[,i]]<-beta
          gamm<-stats::rgamma(n=1,shape=n/2+1,rate=(S[1,1]+lambdaii)/2)
          Omega[i,i]<-gamm+(t(beta) %*% Omega11inv %*% beta)
          
          # Replacing sigma entries
          OmegaInvTemp<-Omega11inv %*% beta
          Sigma[perms[,i],perms[,i]]<-Omega11inv+(OmegaInvTemp %*% t(OmegaInvTemp))/gamm
          Sigma[perms[,i],i]<-(-OmegaInvTemp/gamm)
          Sigma[i,perms[,i]]<-(-OmegaInvTemp/gamm)
          Sigma[i,i]<-1/gamm
        }
        if(iter %% 100==0 & verbose)
        {
          cat("Total iterations= ",iter, "Iterations since burn in= ", 
              ifelse(iter-burnIn>0,iter-burnIn,0), "\n")
        }
        
        # Save Sigma and Omega:
        SigmaMatList[[iter]]<-Sigma
        OmegaMatList[[iter]]<-Omega
      }
    },
    error=function(e)
    {
      # Save Sigma and Omega:
      SigmaMatList[[iter]]<-NULL
      OmegaMatList[[iter]]<-NULL
    }
  )
  bglObj<-list(Sigmas=SigmaMatList,Omegas=OmegaMatList,lambdas=NULL,burnIn=burnIn)
  class(bglObj)<-"BayesianGLasso"
  return(bglObj)
}