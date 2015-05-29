dirichlet.process.prior.mcmc <- function(y,P0,priors,tune,start,n.mcmc){
  
  library(MCMCpack)  # for Dirichlet distribution functions
  
  a0 <- start$a0
  N <- priors$N  # maximum number of clusters
  q.P0 <- length(P0)
  F <- 1/q.P0
  
#   browser()
  tab <- table(y)
  idx <- as.numeric(names(tmp))
  y.tab <- numeric(q.P0)
  y.tab[idx] <- tmp
  p <- rep(a0*F,q.P0) + y.tab
  
  

  P.save <- matrix(0,q.P0,n.mcmc)
  a0.save <- numeric(n.mcmc)
  k.save <- matrix(0,q.P0,n.mcmc)

  keep <- list(a0=0)

  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    
    ###
    ### Sample a0
    ###
    
#     a0.star <- rnorm(1,a0,tune$a0)
#     if (a0.star >= 0) {
#       p.star <- rep(a0.star*F,q.P0) + y.tab
#       mh.star.a0 <- log(ddirichlet(y.tab/100,p.star))+dgamma(a0.star,priors$a,priors$b,log=TRUE)
#       mh.0.a0 <- log(ddirichlet(y.tab/100,p))+dgamma(a0,priors$a,priors$b,log=TRUE)
#       if (exp(mh.star.a0-mh.0.a0)>runif(1)) {
#         a0 <- a0.star
#         p <- p.star
#         keep$a0 <- keep$a0+1
#       }  
#     }
    
    ###
    ### Sample k (cluster assignments for observations)
    ###
    
    p 
    
    
    ###
    ### Sample v
    ###

    v <- c(rbeta(N-1,1,a0),1)
    pie <- v*c(1,cumprod((1-v[-N])))
    
    
    
    ###
    ### Sample P
    ###
    
    P <- rdirichlet(1,p)
    
      
    ###
    ###  Save samples 
    ###

    P.save[, k] <- P
    a0.save[k] <- a0
  }
  
  ###
  ### Write output
  ###
  
  keep$a0 <- keep$a0/n.mcmc
  list(P=P.save,a0=a0.save,keep=keep,n.mcmc=n.mcmc)
  
}