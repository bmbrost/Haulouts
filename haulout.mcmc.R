haulout.mcmc <- function(s,mu,S.tilde,S,priors=list(alpha=1,beta=1,r=0.1,q=0.1),tune=list(mu=0.1),n.mcmc){
  
  ###
  ###  Setup Variables 
  ###
  
  T <- length(mu)
  z.save <- matrix(0,T,n.mcmc)
  mu.save <- matrix(0,T,n.mcmc)
  p.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)
  
  fS.tilde <- 1/dist(S.tilde) # Density of S.tilde
  fS <- 1/dist(S) # Density of S
    
  ###
  ###  Priors and starting values 
  ###
  
  # Starting values for mu
  mu <- s
  idx <- which(s<S.tilde[1])
  mu[idx] <-  S.tilde[1]
  idx <- which(s>S[2])
  mu[idx] <-  S[2]
  
  #Starting values for p and z
  idx <- mu>S.tilde[1]&mu<S.tilde[2] #mu located within intersection(S,S.tilde)
  p <- sum(idx)/T
  z <- numeric(T)
  z[idx] <- 1
  
  sigma <- 0.5
  
  keep <- list(mu=0,sigma=0)
  
  ###
  ### Begin MCMC look
  ###
  
  for(k in 1:n.mcmc){
    if(k%%100==0) cat(k,"");flush.console()
    
    ###
    ### Sample z (haul-out indicator variable)
    ###
    
    idx <- mu>S.tilde[1]&mu<S.tilde[2] #mu located within intersection(S,S.tilde)
    n.tmp <- sum(idx)
    z[!idx] <- 0
    #    browser()
    p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fS)    
    z[idx] <- rbinom(n.tmp,1,p.tmp)
    
    ###
    ### Sample mu (true location)
    ###
    
    #browser()
    
    # Using Gibbs
    mu.star <- rnorm(T,s,sigma) # 'Proposed' mu aren't necessarily in S
    
    # Update mu[t] for z[t]==1
    idx <- which(z==1&mu.star>S.tilde[1]&mu.star<S.tilde[2]) # z[t]==1 and mu.star \in S.tilde
    mu[idx] <- mu.star[idx]

    # Update mu[t] for z[t]==0
    idx <- which(z==0&mu.star>S.tilde[1]&mu.star<S[2]) # z[t]==0 and mu.star \in S
    mu[idx] <- mu.star[idx]
    
    ###
    ### Sample sigma (observation error)
    ###
  
    # Use rigamma {pscl}  
    r.tmp <- sum((s-mu)^2)/2+priors$r
    q.tmp <- T/2+priors$q
    s2 <- rigamma(1,q.tmp,r.tmp)
    sigma <- sqrt(s2)

    # Use rgamma {base}
    #r.tmp <- 1/(sum((s-mu)^2)/2+1/priors$r)
    #q.tmp <- T/2+priors$q
    #s2 <- 1/rgamma(1,q.tmp,,r.tmp)
    #sigma <- sqrt(s2)

    # Use M-H update
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>0&sigma.star<2){
      mh.star.sigma <- sum(dnorm(s,mu,sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(s,mu,sigma,log=TRUE))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
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
    mu.save[,k] <- mu
  }
    
  ###
  ###  Write Output 
  ###

  keep$mu <- keep$mu/(n.mcmc*T)
  keep$sigma <- keep$sigma/n.mcmc
  cat(paste("\nmu acceptance rate:",keep$mu))  
  cat(paste("\nsigma acceptance rate:",keep$sigma))    
  list(z=z.save,p=p.save,mu=mu.save,sigma=sigma.save,keep=keep,n.mcmc=n.mcmc)
}
