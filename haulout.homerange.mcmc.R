haulout.homerange.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc){
  
  library(pscl) #Inverse-gamma random number generator 
  library(mvtnorm)
  library(msm)  #Truncated normal density (dtnorm)
  
  ###
  ###  Setup Variables 
  ###
  
  T <- nrow(mu)
  z.save <- matrix(0,T,n.mcmc)
  mu.save <- array(0,dim=c(T,2,n.mcmc))
  mu.0.save <- matrix(0,n.mcmc,2)
  p.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)
  sigma.mu.save <- numeric(n.mcmc)
  
  # browser()    
  fS.tilde <- 1/(max(dist(S.tilde[,1]))*max(dist(S.tilde[,2]))) # Density of S.tilde
  fS <- 1/(max(dist(S[,1]))*max(dist(S[,2]))) # Density of S

  ###
  ###  Priors and starting values 
  ###
  
  # Starting values for mu
  mu <- start$mu
  sigma <- start$sigma
  mu.0 <- start$mu.0
  sigma.mu <- start$sigma.mu

#   idx <- which(s<S.tilde[1])
#   mu[idx] <-  S.tilde[1]
#   idx <- which(s>S[2])
#   mu[idx] <-  S[2]
  

  #Starting values for p and z
  idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  p <- sum(idx)/T
  z <- numeric(T)
  z[idx] <- 1
  
  keep <- list(mu=0,sigma=0,mu.0=0,sigma.mu=0)
  
  ###
  ### Begin MCMC loop
  ###
  
  for(k in 1:n.mcmc){
    if(k%%100==0) cat(k,"");flush.console()

    ###
    ### Sample z (haul-out indicator variable)
    ###
#     browser()
    idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
      mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    n.tmp <- sum(idx)
    z[!idx] <- 0
    fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
      dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    z[idx] <- rbinom(n.tmp,1,p.tmp)
    
    ###
    ### Sample mu (true location)
    ###
    
#     browser()
    
    # Using full-conditional distribution for mu
    mu.star <- matrix(rnorm(T*2,s,sigma),T,2,byrow=FALSE) # 'Proposed' mu aren't necessarily in S

    # Update mu[t] for z[t]==1
    idx <- which(z==1&mu.star[,1]>S.tilde[1,1]&mu.star[,1]<S.tilde[2,1]&
      mu.star[,2]>S.tilde[1,2]&mu.star[,2]<S.tilde[3,2]) # z[t]==1 and mu.star in S.tilde
    mu[idx,] <- mu.star[idx,]

    # Update mu[t] for z[t]==0
    b <- s%*%solve(sigma^2*diag(2))+matrix(mu.0%*%solve(sigma.mu^2*diag(2)),T,2,byrow=TRUE)
    A <- solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2))    
    A.inv <- solve(A)
    mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))
    mu.star <- cbind(rnorm(T,mu.tmp[,1],sqrt(A.inv[1,1])),rnorm(T,mu.tmp[,2],sqrt(A.inv[2,2])))
    idx <- which(z==0&mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2]) # z[t]==0 and mu.star in S
    mu[idx,] <- mu.star[idx,]

  
    ###
    ### Sample sigma (observation error)
    ###
  
    # Use rigamma {pscl}  
#     r.tmp <- sum((s-mu)^2)/2+priors$r
#     q.tmp <- T+priors$q
#     s2 <- rigamma(1,q.tmp,r.tmp)
#     sigma <- sqrt(s2)

    # Use rgamma {base}
    #r.tmp <- 1/(sum((s-mu)^2)/2+1/priors$r)
    #q.tmp <- T/2+priors$q
    #s2 <- 1/rgamma(1,q.tmp,,r.tmp)
    #sigma <- sqrt(s2)

    # Use M-H update
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$a&sigma.star<priors$b){
      mh.star.sigma <- sum(dnorm(s,mu,sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(s,mu,sigma,log=TRUE))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
      } 
    }

    ###
    ### Sample mu.0 (homerange center)
    ###

	# Could sample using Gibbs and then throw out proposals not in S.tilde?

    mu.0.star <- rnorm(2,mu.0,tune$mu.0)
    if(mu.0.star[1]>S.tilde[1,1]&mu.0.star[1]<S.tilde[2,1]& #Reject proposals for mu.0 not in S.tilde
      mu.0.star[2]>S.tilde[1,2]&mu.0.star[2]<S.tilde[3,2]){
      idx <- which(z==0) #Update using z==0 only
      mh.star.mu.0 <- sum(dtnorm(mu[idx,1],mu.0.star[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
        dtnorm(mu[idx,2],mu.0.star[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      mh.0.mu.0 <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
        dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      if(exp(mh.star.mu.0-mh.0.mu.0)>runif(1)){
        mu.0 <- mu.0.star
        keep$mu.0 <- keep$mu.0+1
      } 
    }    

    
    ###
    ### Sample sigma.mu (disperson around homerange center)
    ###

    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    if(sigma.mu.star>0){
      mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu.star,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
        dtnorm(mu[idx,2],mu.0[2],sigma.mu.star,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
        dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        sigma.mu <- sigma.mu.star
        keep$sigma.mu <- keep$sigma.mu+1
      } 
    }    
    
    ###
    ### Sample p (probability of hauled out)
    ###
    
    sumz <- sum(z)
    p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)

    ###
    ###  Save samples 
    ###
    
    p.save[k] <- p
    sigma.save[k] <- sigma
    z.save[,k] <- z
    mu.save[,,k] <- mu
    mu.0.save[k,] <- mu.0
    sigma.mu.save[k] <- sigma.mu
  }
    
  ###
  ###  Write Output 
  ###

  keep$mu <- keep$mu/(n.mcmc*T)
  keep$sigma <- keep$sigma/n.mcmc
  keep$mu.0 <- keep$mu.0/n.mcmc
  keep$sigma.mu <- keep$sigma.mu/n.mcmc
  cat(paste("\nmu acceptance rate:",keep$mu))  
  cat(paste("\nsigma acceptance rate:",keep$sigma))    
  cat(paste("\nmu.0 acceptance rate:",keep$mu.0))    
  cat(paste("\nsigma.mu acceptance rate:",keep$sigma.mu))    
  list(z=z.save,p=p.save,mu=mu.save,sigma=sigma.save,mu.0=mu.0.save,sigma.mu=sigma.mu.save,keep=keep,n.mcmc=n.mcmc)
}
