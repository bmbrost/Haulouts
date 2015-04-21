haulout.mcmc <- function(mu,S.tilde,S,priors=list(alpha=1,beta=1),n.mcmc){
  
  ###
  ###  Setup Variables 
  ###
  
  T <- length(mu)
  
  z.save <- matrix(0,T,n.mcmc)
  p.save <- numeric(n.mcmc)
  
  r1 <- 1/dist(S.tilde)
  r2 <- 1/dist(S)
  
  ###
  ###  Priors and starting values 
  ###
  
  p <- 0.5
  z <- rbinom(T,1,p)
  
  ###
  ### Begin MCMC look
  ###
  
  for(k in 1:n.mcmc){
    if(k%%100==0) cat(k,"");flush.console()
    
    ###
    ### Sample p
    ###
    
    sumz <- sum(z)
    p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)
    
    
    ###
    ### Sample z
    ###
  
    p.tmp <- (p*r1)/(p*r1+(1-p)*r2)
    p.tmp <- p/(p+(1-p))
    z <- rbinom(T,1,p.tmp)
#     z <- numeric(T)
#     z[mu>H[1]&mu<H[2]] <- 1
    
    ###
    ###  Save samples 
    ###
    
    p.save[k] <- p
    z.save[,k] <- z
  }
    
  ###
  ###  Write Output 
  ###
  
  list(z=z.save,p=p.save)
}
