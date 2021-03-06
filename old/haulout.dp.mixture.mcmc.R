haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  t.start <- Sys.time()

  ###
  ### Libraries and Subroutines
  ###
  
  library(MCMCpack)  # for Dirichlet distribution functions
  library(data.table)  # for tabulating and summing
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  library(doParallel)  # for parallel processing
  library(foreach)  # for parallel processing  
  # library(mvtnorm)  # for multivariate normal density
  library(msm)  # for truncated normal density
  library(data.table)  # for data.table tabulation functionality
  
	get.mu.0 <- function(x,h.idx,z,s,Sigma.inv,Sigma.mu.inv){
		# browser()
		idx.0 <- which(h.idx==x&z==0)
		idx.1 <- which(h.idx==x&z==1)
		n.0 <- length(idx.0)
		n.1 <- length(idx.1)
		b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(s[idx.0,]%*%(Sigma.inv+Sigma.mu.inv))
		A.inv <- solve(n.1*Sigma.inv+n.0*(Sigma.inv+Sigma.mu.inv))
		rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
	}

  
  ###
  ###  Create cluster for parallel processing
  ###
  
  #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  #   mcoptions <- list(preschedule=TRUE)
  #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
 # browser() 
  ###
  ### Starting values and priors
  ###
  
  a0 <- start$a0
  z <- start$z
  sigma <- start$sigma
  sigma.mu <- start$sigma.mu
  pie <- start$pie
  H <- priors$H
  p <- start$p

  Sigma.inv <- solve(sigma^2*diag(2))
  Sigma.mu.inv <- solve(sigma.mu^2*diag(2))

 #Starting values for p and z
  # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  # p <- sum(idx)/T
  # z <- numeric(T)
  # z[idx] <- 1

  
  ###
  ###  Setup Variables 
  ###
  
  # browser() 
	T <- nrow(s)  # number of observations
	mu.0 <- unique(h)  # unique cluster location s
	n.cls <- nrow(mu.0)  # number of clusters
	h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator

	# Tabulate cluster membership with data.table and setdiff	
	# dt.h.idx <- as.data.table(h.idx)
	# dt.tab.cls <- dt.h.idx[,.N,by=h.idx]
	# setkey(dt.tab.cls,N)
	# idx.cls <- rev(dt.tab.cls[,h.idx])
	# samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order in which clusters are sampled

	# Tabulate cluster membership with data.table but not setdiff
	# h.idx <- c(h.idx,1:H)
	# dt.h.idx <- as.data.table(h.idx)
	# dt.tab.cls <- dt.h.idx[,.N,by=h.idx]
	# dt.tab.cls[,N:=N-1]
	# setkey(dt.tab.cls,N)
	# idx.cls <- rev(dt.tab.cls[N>0,h.idx])
	# samp.cls <- rev(dt.tab.cls[,h.idx])  # order in which clusters are sampled

	 
	# Propose values for mu.0
	n.cls.star <- H-n.cls  # number of new clusters to propose
	mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star


  ###
  ### Create receptacles for output
  ###
  
  h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  a0.save <- numeric(n.mcmc)  # concentration parameter
  sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  p.save <- numeric(n.mcmc)  # probability of being hauled-out
  z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  n.cls.save <- numeric(n.mcmc)
  
  keep <- list(sigma=0,sigma.mu=0)

    
  ###
  ### Begin MCMC loop
  ###
  
  for (k in 1:n.mcmc) {
    if(k%%1000==0) cat(k,"");flush.console()

	###
	### Dirichlet process parameters
	###

	# Update follows the blocked Gibbs sampler of Ishwaran and James (2001) and
    # Gelman et al. (2014), Section 23.3

	# Note: sampling order matters here. Cluster parameters must be sampled in same 
	# order as pie, i.e., sorted by decreasing membership

	   	### Sample h.t (cluster assignment indicator)
# browser()
		# Sampled with truncated normal density
	   	h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
			prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
			dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
			(dnorm(s[x,1],mu.0[samp.cls,1],sqrt(sigma^2+sigma.mu^2))*
			dnorm(s[x,2],mu.0[samp.cls,2],sqrt(sigma^2+sigma.mu^2)))^(1-z[x])))

		# Tabulate cluster membership with data.table but not setdiff
		# h.idx <- c(h.idx,1:H)
		# dt.h.idx <- as.data.table(h.idx)
		# dt.tab.cls <- dt.h.idx[,.N,by=h.idx]
		# dt.tab.cls[,N:=N-1]
		# setkey(dt.tab.cls,N)
		# n.cls <- dt.tab.cls[N>0,.N]
		# idx.cls <- rev(dt.tab.cls[N>0,h.idx])
		# samp.cls <- rev(dt.tab.cls[,h.idx])  # order of decreasing membership
	
		# Tabulate cluster membership with data.table and setdiff	
		# dt.h.idx <- as.data.table(h.idx)
		# dt.tab.cls <- dt.h.idx[,.N,by=h.idx]
		# setkey(dt.tab.cls,N)
		# n.cls <- dt.tab.cls[,.N]
		# idx.cls <- rev(dt.tab.cls[,h.idx])
		# samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership

		# Tabulate cluster membership with base functions
		tab.cls <- table(h.idx)  # tabulate cluster membership
		n.cls <- length(tab.cls)  # number of clusters
		ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
		tab.cls <- tab.cls[ord]
		idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
		samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership
  
 	    ### Sample pie (stick-breaking process)
    
	    # Use for data.table functionality
		# tab.cls.tmp <- c(rev(dt.tab.cls[,N]),rep(0,H-n.cls-1)) 
		
		# Use for base functionality
 	  	tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))  # membership in decreasing order
	    v <- c(rbeta(H-1,1+tab.cls.tmp,a0+T-cumsum(tab.cls.tmp)),1)  # stick-breaking weights
	    pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	    ### Sample a0 (concentration parameter); See Gelman section 23.3
       
    	a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
# a0 <- start$a0

	###
    ### Sample mu.0 (true location of occupied clusters)
    ###
	   
	# Sampling order does not matter here
	# browser()	

	# Use for data.table functionality
	# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0
		# get.mu.0(x,dt.h.idx[1:T,h.idx],z,s,mu,sigma,sigma.mu,S.tilde)))
	
	# Use for base functionality	
	mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0	
		get.mu.0(x,h.idx,z,s,Sigma.inv,Sigma.mu.inv)))  
	idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	
	# if(length(idx)>0){
		mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # accept proposals in S.tilde
	# }
	n.cls.star <- H-n.cls  # number of new clusters to propose
	mu.0[-idx.cls,] <- cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      runif(n.cls.star,S.tilde[1,2],S.tilde[3,2]))  # update mu.0 with mu.star


    ###
    ### Sample sigma (observation error)
    ###
# browser()
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
		idx <- z==1

		# mu.0.tmp1 <- mu.0[h.idx[idx],] 
		# mu.0.tmp0 <- mu.0[h.idx[-idx],] 

		mu.0.tmp <- mu.0[h.idx,]
    	mh.star.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma.star,log=TRUE)+
	    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma.star,log=TRUE))#+
	    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sqrt(sigma.star^2+sigma.mu^2),log=TRUE)+
	    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sqrt(sigma.star^2+sigma.mu^2),log=TRUE))
	    mh.0.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma,log=TRUE)+
	    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma,log=TRUE))#+
	    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sqrt(sigma^2+sigma.mu^2),log=TRUE)+
	    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sqrt(sigma^2+sigma.mu^2),log=TRUE))
    	# mh.star.sigma <- sum(sapply(1:T, function(x) 
      		# dmvnorm(s[x,],mu[x,],sigma.star^2*diag(2),log=TRUE)))
      	# mh.0.sigma <- sum(sapply(1:T, function(x) 
      		# dmvnorm(s[x,],mu[x,],sigma^2*diag(2),log=TRUE)))
	    if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
    	    sigma <- sigma.star
        	Sigma.inv <- solve(sigma^2*diag(2))
	        keep$sigma <- keep$sigma+1
    	} 
    }
# sigma <- start$sigma


    ###
    ### Sample sigma.mu (disperson around homerange center)
    ###
# browser()

    # Sample with truncated normal density
    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  idx <- which(z==0)
	  mu.0.tmp <- mu.0[h.idx[idx],]
      # print(idx)
      mh.star.sigma.mu <- 
      	  sum(dnorm(s[idx,1],mu.0.tmp[,1],sqrt(sigma^2+sigma.mu.star^2),log=TRUE)+
	      dnorm(s[idx,2],mu.0.tmp[,2],sqrt(sigma^2+sigma.mu.star^2),log=TRUE))
      mh.0.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sqrt(sigma^2+sigma.mu^2),log=TRUE)+
	      dnorm(s[idx,2],mu.0.tmp[,2],sqrt(sigma^2+sigma.mu^2),log=TRUE))
      if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        sigma.mu <- sigma.mu.star
        Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
        keep$sigma.mu <- keep$sigma.mu+1
      } 
    }
# sigma.mu <- start$sigma.mu
	

	###
    ### Sample z (haul-out indicator variable)
    ###

	p.tmp1 <- p*dnorm(s[,1],mu.0[h.idx,1],sigma,log=FALSE)*
		dnorm(s[,2],mu.0[h.idx,2],sigma,log=FALSE)
	p.tmp2 <- (1-p)*dnorm(s[,1],mu.0[h.idx,1],sqrt(sigma^2+sigma.mu^2),log=FALSE)*
		dnorm(s[,2],mu.0[h.idx,2],sqrt(sigma^2+sigma.mu^2),log=FALSE)
	p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
	z <- rbinom(T,1,p.tmp)
# z <- start$z

	
	###
    ### Sample p (probability of hauled out)
    ###
    
    sumz <- sum(z)
    p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)

	
    ###
    ###  Save samples 
    ###

	# Use with data.table functionality
	# h.idx.save[,k] <- dt.h.idx[1:T,h.idx]
	# h.save[,,k] <- mu.0[dt.h.idx[1:T,h.idx],]
	
	# Use with base functionality
	h.idx.save[,k] <- h.idx 
	h.save[,,k] <- mu.0[h.idx,]
	mu.0.save[,,k] <- mu.0
    a0.save[k] <- a0    
    sigma.save[k] <- sigma
    sigma.mu.save[k] <- sigma.mu
    p.save[k] <- p
	z.save[,k] <- z
	n.cls.save[k] <- n.cls
  }
  
  ###
  ### Write output
  ###
  
  keep$sigma <- keep$sigma/n.mcmc
  keep$sigma.mu <- keep$sigma.mu/n.mcmc
  cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
  cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
  cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,a0=a0.save,sigma=sigma.save,
  	sigma.mu=sigma.mu.save,p=p.save,z=z.save,n.cls=n.cls.save,keep=keep,n.mcmc=n.mcmc)
  
}