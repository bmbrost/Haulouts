haulout.2d.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc){
  
  library(pscl) #Inverse-gamma random number generator 
  
  ###
  ###  Setup Variables 
  ###
  
  T <- nrow(mu)
  z.save <- matrix(0,T,n.mcmc)
  mu.save <- array(0,dim=c(T,2,n.mcmc))
  p.save <- numeric(n.mcmc)
  sigma.save <- numeric(n.mcmc)
# browser()    
  fS.tilde <- 1/(max(dist(S.tilde[,1]))*max(dist(S.tilde[,2]))) # Density of S.tilde
  fS <- 1/(max(dist(S[,1]))*max(dist(S[,2]))) # Density of S

  ###
  ###  Priors and starting values 
  ###
  
  # Starting values for mu
  mu <- start$mu
#   idx <- which(s<S.tilde[1])
#   mu[idx] <-  S.tilde[1]
#   idx <- which(s>S[2])
#   mu[idx] <-  S[2]
  

  #Starting values for p and z
  idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  p <- sum(idx)/T
  z <- numeric(T)
  z[idx] <- 1
  
  sigma <- start$sigma
  keep <- list(mu=0,sigma=0)
  
  ###
  ### Begin MCMC look
  ###
  
  for(k in 1:n.mcmc){
    if(k%%100==0) cat(k,"");flush.console()
    
    ###
    ### Sample z (haul-out indicator variable)
    ###
    
    idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
      mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    n.tmp <- sum(idx)
    z[!idx] <- 0
    p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fS)    
    z[idx] <- rbinom(n.tmp,1,p.tmp)
    
    ###
    ### Sample mu (true location)
    ###
    
    # browser()
    
    ### Using full-conditional distribution for mu
    mu.star <- matrix(rnorm(T*2,s,sigma),T,2,byrow=FALSE) # 'Proposed' mu aren't necessarily in S
      
    # Update mu[t] for z[t]==1 using full-conditional distribution for mu
    idx <- which(z==1&mu.star[,1]>S.tilde[1,1]&mu.star[,1]<S.tilde[2,1]&
      mu.star[,2]>S.tilde[1,2]&mu.star[,2]<S.tilde[3,2]) # z[t]==1 and mu.star in S.tilde
    mu[idx,] <- mu.star[idx,] 
    
    # Update mu[t] for z[t]==0 using full-conditional distribution for mu
    idx <- which(z==0&mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
       mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2]) # z[t]==0 and mu.star in S
    mu[idx,] <- mu.star[idx,]

    
    ### Using MH update
#     mu.star <- matrix(rnorm(T*2,mu,tune$mu),T,2,byrow=FALSE)
#     
#     # Update mu[t] for z[t]==1 using full-conditional distribution for mu
#     idx <- which(z==1&mu.star[,1]>S.tilde[1,1]&mu.star[,1]<S.tilde[2,1]&
#       mu.star[,2]>S.tilde[1,2]&mu.star[,2]<S.tilde[3,2]) # z[t]==1 and mu.star in S.tilde
#     mh.star.mu <- rowSums(dnorm(s[idx,],mu.star[idx,],sigma,log=TRUE))
#     mh.0.mu <- rowSums(dnorm(s[idx,],mu[idx,],sigma,log=TRUE))
#     idx.keep <- idx[exp(mh.star.mu-mh.0.mu)>runif(length(idx))]
#     mu[idx.keep,] <- mu.star[idx.keep,]
#     keep$mu <- keep$mu+length(idx.keep)
#     
#     # Update mu[t] for z[t]==0 using full-conditional distribution for mu
#     idx <- which(z==0&mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
#       mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2]) # z[t]==0 and mu.star in S
#     mh.star.mu <- rowSums(dnorm(s[idx,],mu.star[idx,],sigma,log=TRUE))
#     mh.0.mu <- rowSums(dnorm(s[idx,],mu[idx,],sigma,log=TRUE))
#     idx.keep <- idx[exp(mh.star.mu-mh.0.mu)>runif(length(idx))]
#     mu[idx.keep,] <- mu.star[idx.keep,]
#     keep$mu <- keep$mu+length(idx.keep)
#   

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
