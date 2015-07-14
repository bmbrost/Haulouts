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
  
	get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		idx.0 <- which(h.idx==x&z==0)
		idx.1 <- which(h.idx==x&z==1)
		n.0 <- length(idx.0)
		n.1 <- length(idx.1)
		Sigma.inv <- solve(sigma^2*diag(2))
		Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(mu[idx.0,]%*%Sigma.mu.inv)
		A <- n.1*Sigma.inv+n.0*Sigma.mu.inv
		A.inv <- solve(A)
		mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		mu.0.tmp
	}
  
  
  ###
  ###  Create cluster for parallel processing
  ###
  
  #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  #   mcoptions <- list(preschedule=TRUE)
  #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  ###
  ### Starting values and priors
  ###
  
  a0 <- start$a0
  z <- start$z
  sigma <- start$sigma
  sigma.mu <- start$sigma.mu
  pie <- start$pie
  H <- priors$H
  mu <- start$mu
  p <- start$p

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

	# Tabulate cluster membership with base functions
	tab.cls <- table(h.idx)  # tabulate cluster membership
	ord <- order(tab.cls,decreasing=TRUE) # order of clusters by membership
	tab.cls <- tab.cls[ord]  # ordered largest to smallest
	idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters in order
	samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order in which clusters are sampled
 
	# Propose values for mu.0
	n.cls.star <- H-n.cls  # number of new clusters to propose
	mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star


  ###
  ### Create receptacles for output
  ###
  
  h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  mu.save <- array(0,dim=c(T,2,n.mcmc))  # true animal locations
  a0.save <- numeric(n.mcmc)  # concentration parameter
  sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  p.save <- numeric(n.mcmc)  # probability of being hauled-out
  z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  
  keep <- list(mu=0,sigma=0,sigma.mu=0)

    
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
			(dtnorm(rep(mu[x,1],H),mu.0[samp.cls,1],sigma.mu,
			lower=min(S[,1]),upper=max(S[,1]))*
			dtnorm(rep(mu[x,2],H),mu.0[samp.cls,2],sigma.mu,lower=min(S[,2]),
			upper=max(S[,2])))^(1-z[x])))

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
a0 <- start$a0


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
		get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  
	idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	
	mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # accept proposals in S.tilde
	n.cls.star <- H-n.cls  # number of new clusters to propose
	mu.0[-idx.cls,] <- cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      runif(n.cls.star,S.tilde[1,2],S.tilde[3,2]))  # update mu.0 with mu.star


    ###
    ### Sample mu (true location of individual)
    ###
    
# browser()
	# Sampling order does not matter here
	
    # Update mu[t] for z[t]==1, locations at haul-out sites
	idx <- which(z==1)  # locations at haul-out
	h.idx.tmp <- h.idx[idx]
	mu[idx,] <- mu.0[h.idx.tmp,]  # mu is haul-out site for hauled-out individuals

    # Update mu[t] for z[t]==0, at-sea locations
	idx <- which(z==0)  # at-sea locations
	h.idx.tmp <- h.idx[idx]
    b <- s[idx,]%*%solve(sigma^2*diag(2))+mu.0[h.idx.tmp,]%*%solve(sigma.mu^2*diag(2))
    A.inv <- solve(solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2)))  # var-cov matrix    
    mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))  # mean matrix
	T.0 <- nrow(mu.tmp)   
    mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))  # proposals for mu
    idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	mu[idx[idx.tmp],] <- mu.star[idx.tmp,]
# mu <- start$mu
    
    ###
    ### Sample sigma (observation error)
    ###
# browser()
    sigma.star <- rnorm(1,sigma,tune$sigma)
    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      mh.star.sigma <- sum(dnorm(s[,1],mu[,1],sigma.star,log=TRUE)+
	      dnorm(s[,2],mu[,2],sigma.star,log=TRUE))
      mh.0.sigma <- sum(dnorm(s[,1],mu[,1],sigma,log=TRUE)+
	      dnorm(s[,2],mu[,2],sigma,log=TRUE))
      # mh.star.sigma <- sum(sapply(1:T, function(x) 
      	# dmvnorm(s[x,],mu[x,],sigma.star^2*diag(2),log=TRUE)))
      # mh.0.sigma <- sum(sapply(1:T, function(x) 
      	# dmvnorm(s[x,],mu[x,],sigma^2*diag(2),log=TRUE)))
      if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        sigma <- sigma.star
        keep$sigma <- keep$sigma+1
      } 
    }
sigma <- start$sigma
    
    ###
    ### Sample sigma.mu (disperson around homerange center)
    ###
# browser()

    # Sample with truncated normal density
    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  idx <- which(z==0)
	  mu.0.tmp <- mu.0[h.idx[idx],]
      mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,
	      lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,
	      lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,
	      lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,
	      lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        sigma.mu <- sigma.mu.star
        keep$sigma.mu <- keep$sigma.mu+1
      } 
    }
sigma.mu <- start$sigma.mu

	###
    ### Sample z (haul-out indicator variable)
    ###

# browser()
    idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
	n.tmp <- sum(idx)
	z[!idx] <- 0
	h.idx.tmp <- h.idx[idx]
	# mu.0.tmp <- mu.0[h.idx[idx],]
	
# points(mu[idx,],cex=p.tmp1)
# points(s[idx,],cex=p.tmp1,col=2)
# plot(mu,col=idx+1)

# browser()

	# p.tmp1 <- p*dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		# dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE)

	p.tmp1 <- p*dnorm(s[idx,1],mu.0[h.idx.tmp,1],sigma,log=FALSE)*
		dnorm(s[idx,2],mu.0[h.idx.tmp,2],sigma,log=FALSE)
	# p.tmp1 <- p

# plot(p.tmp1,col=z[idx]+1)
# head(h.idx)
# hist(p.tmp,breaks=50)
# points(mu[idx,], cex=p.tmp)
# points(s[idx,], cex=p.tmp,col=2)

	p.tmp2 <- (1-p)*dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE)*
		dtnorm(mu[idx,1],mu.0[h.idx.tmp,1],sigma.mu,
			lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
		dtnorm(mu[idx,2],mu.0[h.idx.tmp,2],sigma.mu,
			lower=min(S[,2]),upper=max(S[,2]),log=FALSE)

	# p.tmp2 <- (1-p)*(0.5*(dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		# dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE))+
		# (0.5*dtnorm(mu[idx,1],mu.0[h.idx.tmp,1],sigma.mu,
			# lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
		# dtnorm(mu[idx,2],mu.0[h.idx.tmp,2],sigma.mu,
			# lower=min(S[,2]),upper=max(S[,2]),log=FALSE)))

	# p.tmp2 <- (1-p)*dnorm(s[idx,1],mu.0[h.idx.tmp,1],sqrt(sigma^2+sigma.mu^2),log=FALSE)*
		# dnorm(s[idx,2],mu.0[h.idx.tmp,2],sqrt(sigma^2+sigma.mu^2),log=FALSE)

# plot(p.tmp20,p.tmp2,col=z[idx]+1)
# plot(mu[idx,1],mu.0[h.idx.tmp,1],col=z[idx]+1)
# abline(a=0,b=1)
# plot(p.tmp,col=z[idx]+1)
# plot(p.tmp,p.tmp2,col=z[idx]+1)
# plot(mu[idx,],col=idx+1,cex=p.tmp+1)
# points(mu.0[h.idx.tmp,],col=3,pch=19,cex=0.5)

	p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
	z[idx] <- rbinom(n.tmp,1,p.tmp)

# z <- start$z


	# ###
    # ### Sample p (probability of hauled out)
    # ###
    
    sumz <- sum(z)
    p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)
# p <- start$p

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
    mu.save[,,k] <- mu
    p.save[k] <- p
	z.save[,k] <- z
  }
  
  ###
  ### Write output
  ###
  
  keep$sigma <- keep$sigma/n.mcmc
  keep$sigma.mu <- keep$sigma.mu/n.mcmc
  cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  cat(paste("\nsigma.mu acceptance rate:",keep$sigma.mu)) 
  cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,mu=mu.save,a0=a0.save,sigma=sigma.save,
  	sigma.mu=sigma.mu.save,p=p.save,z=z.save,keep=keep,n.mcmc=n.mcmc)
  
}





# haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  # t.start <- Sys.time()

  # ###
  # ### Libraries and Subroutines
  # ###
  
  # library(MCMCpack)  # for Dirichlet distribution functions
  # library(data.table)  # for tabulating and summing
  # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # library(doParallel)  # for parallel processing
  # library(foreach)  # for parallel processing  
  # # library(mvtnorm)  # for multivariate normal density
  # library(msm)  # for truncated normal density
  
	# # get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# # idx.0 <- which(h.idx==x&z==0)
		# # idx.1 <- which(h.idx==x&z==1)
		# # n.0 <- length(idx.0)
		# # n.1 <- length(idx.1)
		# # b <- colSums(s[idx.1,]%*%solve(sigma^2*diag(2)))+
			# # colSums(mu[idx.0,]%*%solve(sigma.mu^2*diag(2)))
		# # A <- n.1*solve(sigma^2*diag(2))+n.0*solve(sigma.mu^2*diag(2))
		# # A.inv <- solve(A)
		# # # mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# # mu.0.tmp
	# # }
  
  
  # ###
  # ###  Create cluster for parallel processing
  # ###
  
  # #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  # #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  # #   mcoptions <- list(preschedule=TRUE)
  # #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # ###
  # ### Starting values and priors
  # ###
  
  # a0 <- start$a0
  # z <- start$z
  # sigma <- start$sigma
  # sigma.mu <- start$sigma.mu
  # pie <- start$pie
  # H <- priors$H
  # mu <- start$mu
  # p <- start$p

 # #Starting values for p and z
  # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  # # p <- sum(idx)/T
  # # z <- numeric(T)
  # # z[idx] <- 1

  
  # ###
  # ###  Setup Variables 
  # ###
  
  # # browser() 
  # T <- nrow(s)  # number of observations
  # mu.0 <- unique(h)  # unique cluster locations
# # browser()
  # n.cls <- nrow(mu.0)  # number of clusters
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  # tab.cls <- table(h.idx)  # tabulate cluster membership
 
  # # Order by decreasing membership
  # ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
  # tab.cls <- tab.cls[ord]
  # idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
 
  # # Propose values for mu.0
  # n.cls.star <- H-n.cls  # number of new clusters to propose
  # mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star
  # samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership  

  # ###
  # ### Create receptacles for output
  # ###
  
  # h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
# h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  # mu.save <- array(0,dim=c(T,2,n.mcmc))  # true animal locations
  # a0.save <- numeric(n.mcmc)  # concentration parameter
  # sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  # sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  # p.save <- numeric(n.mcmc)  # probability of being hauled-out
  # z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  # mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  
  # keep <- list(mu=0,sigma=0,sigma.mu=0)

    
  # ###
  # ### Begin MCMC loop
  # ###
  
  # for (k in 1:n.mcmc) {
    # if(k%%1000==0) cat(k,"");flush.console()

	# ###
	# ### Dirichlet process parameters
	# ###

	# # Update follows the blocked Gibbs sampler of Ishwaran and James (2001) and
    # # Gelman et al. (2014), Section 23.3

	# # Note: sampling order matters here. Clusters must be sampled in same 
	# # order as pie, i.e., sorted by decreasing membership

	   	# ### Sample h.t (cluster assignment indicator)
# # browser()
	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # # (dnorm(mu[x,1],mu.0[samp.cls,1],sigma.mu)*
		  # # dnorm(mu[x,2],mu.0[samp.cls,2],sigma.mu))^(1-z[x])))

		# # Sampled with truncated normal density
	
	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]))

	   	# h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # (dtnorm(rep(mu[x,1],H),mu.0[samp.cls,1],sigma.mu,
		  # lower=min(S[,1]),upper=max(S[,1]))*
		  # dtnorm(rep(mu[x,2],H),mu.0[samp.cls,2],sigma.mu,lower=min(S[,2]),
		  # upper=max(S[,2])))^(1-z[x])))
		# tab.cls <- table(h.idx)  # tabulate cluster membership
		# n.cls <- length(tab.cls)  # number of clusters
	
	  	# # Sort in order of decreasing cluster membership
		# ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
		# tab.cls <- tab.cls[ord]
		# idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
		# samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership
  
 	    # ### Sample pie (stick-breaking process)
    
	  	# tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))  # membership in decreasing order
	    # v <- c(rbeta(H-1,1+tab.cls.tmp,a0+T-cumsum(tab.cls.tmp)),1)  # stick-breaking weights
	    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	    # ### Sample a0 (concentration parameter); See Gelman section 23.3
       
    	# a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
# # a0 <- start$a0

    # ###
    # ### Sample mu.0 (true location of occupied clusters)
    # ###
	   
	# # Sampling order does not matter here

	# get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# idx.0 <- which(h.idx==x&z==0)
		# idx.1 <- which(h.idx==x&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# Sigma.inv <- solve(sigma^2*diag(2))
		# Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		# b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(mu[idx.0,]%*%Sigma.mu.inv)
		# A <- n.1*Sigma.inv+n.0*Sigma.mu.inv
		# A.inv <- solve(A)
		# mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# mu.0.tmp
	# }

# # browser()	
	# mu.0.tmp <- t(sapply(idx.cls,function(x) 
		# get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  # proposals for mu.0	

	# # idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		# # mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	

	# # mu.0[idx,] <- mu.0.tmp[idx,]

	# mu.0[idx.cls,] <- mu.0.tmp
	
	# # mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # update mu.0
	# # mu.0[idx.cls,] <- mu.0.tmp  # update mu.0

	# n.cls.star <- H-n.cls  # number of new clusters to propose
	# mu.0[-idx.cls,] <- cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2]))  # update mu.0 with mu.star


    # ###
    # ### Sample mu (true location of individual)
    # ###
    
# # browser()
	# # Sampling order does not matter here
	
    # # Update mu[t] for z[t]==1, locations at haul-out sites

	# idx <- which(z==1)  # locations at haul-out
	# h.idx.tmp <- h.idx[idx]
	# mu[idx,] <- mu.0[h.idx.tmp,]  # mu is haul-out site for hauled-out individuals

    # # Update mu[t] for z[t]==0, at-sea locations
	# # idx <- which(z==0)  # at-sea locations
	# # h.idx.tmp <- h.idx[idx]
    # # b <- s[idx,]%*%solve(sigma^2*diag(2))+mu.0[h.idx.tmp,]%*%solve(sigma.mu^2*diag(2))
    # # A.inv <- solve(solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2)))  # var-cov matrix    
    # # mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))  # mean matrix
	# # T.0 <- nrow(mu.tmp)   
    # # mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	# # rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))  # proposals for mu
    # # idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      # # mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	# # mu[idx[idx.tmp],] <- mu.star[idx.tmp,]
# # mu <- start$mu
    
    # ###
    # ### Sample sigma (observation error)
    # ###
# # browser()
    # sigma.star <- rnorm(1,sigma,tune$sigma)
    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # mh.star.sigma <- sum(dnorm(s[,1],mu[,1],sigma.star,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma.star,log=TRUE))
      # mh.0.sigma <- sum(dnorm(s[,1],mu[,1],sigma,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma,log=TRUE))
      # # mh.star.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma.star^2*diag(2),log=TRUE)))
      # # mh.0.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma^2*diag(2),log=TRUE)))
      # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # sigma <- sigma.star
        # keep$sigma <- keep$sigma+1
      # } 
    # }
# # sigma <- start$sigma
    
    # ###
    # ### Sample sigma.mu (disperson around homerange center)
    # ###
# # browser()
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,log=TRUE))
      # # mh.0.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

    # # Sample with truncated normal density
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

# # sigma.mu <- start$sigma.mu

	# ###
    # ### Sample z (haul-out indicator variable)
    # ###

# # browser()
    # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	# # mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    # # n.tmp <- sum(idx)
    # # z[!idx] <- 0
	# # p.tmp1 <- (p*dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		# # dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE))^z[idx]
	# # # p.tmp2 <- ((1-p)*dnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu)*
		# # # dnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu))^(1-z[idx])
	# # # Sample with truncated normal density
	# # p.tmp2 <- ((1-p)*dtnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu,
		# # lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
		# # dtnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu,
		# # lower=min(S[,2]),upper=max(S[,2]),log=FALSE))^(1-z[idx])
	# # p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
	# # z[idx] <- rbinom(n.tmp,1,p.tmp)

    # # # fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
	    # # # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    # # # p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    # # # z[idx] <- rbinom(n.tmp,1,p.tmp)

# # z <- start$z



	# # ###
    # # ### Sample p (probability of hauled out)
    # # ###
    
    # # sumz <- sum(z)
    # # p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)
# # p <- start$p

    # ###
    # ###  Save samples 
    # ###

	# h.idx.save[,k] <- h.idx 
# h.save[,,k] <- mu.0[h.idx,]
	# mu.0.save[,,k] <- mu.0
    # a0.save[k] <- a0    
    # sigma.save[k] <- sigma
    # sigma.mu.save[k] <- sigma.mu
    # mu.save[,,k] <- mu
    # p.save[k] <- p
	# z.save[,k] <- z
  # }
  
  # ###
  # ### Write output
  # ###
  
  # keep$sigma <- keep$sigma/n.mcmc
  # keep$sigma.mu <- keep$sigma.mu/n.mcmc
  # cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  # cat(paste("\nsigma.mu acceptance rate:",keep$sigma.mu)) 
  # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  # list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,mu=mu.save,a0=a0.save,sigma=sigma.save,
  	# sigma.mu=sigma.mu.save,p=p.save,z=z.save,keep=keep,n.mcmc=n.mcmc)
  
# }




# haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  # t.start <- Sys.time()

  # ###
  # ### Libraries and Subroutines
  # ###
  
  # library(MCMCpack)  # for Dirichlet distribution functions
  # library(data.table)  # for tabulating and summing
  # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # library(doParallel)  # for parallel processing
  # library(foreach)  # for parallel processing  
  # # library(mvtnorm)  # for multivariate normal density
  # library(msm)  # for truncated normal density
  
	# # get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# # idx.0 <- which(h.idx==x&z==0)
		# # idx.1 <- which(h.idx==x&z==1)
		# # n.0 <- length(idx.0)
		# # n.1 <- length(idx.1)
		# # b <- colSums(s[idx.1,]%*%solve(sigma^2*diag(2)))+
			# # colSums(mu[idx.0,]%*%solve(sigma.mu^2*diag(2)))
		# # A <- n.1*solve(sigma^2*diag(2))+n.0*solve(sigma.mu^2*diag(2))
		# # A.inv <- solve(A)
		# # # mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# # mu.0.tmp
	# # }
  
  
  # ###
  # ###  Create cluster for parallel processing
  # ###
  
  # #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  # #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  # #   mcoptions <- list(preschedule=TRUE)
  # #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # ###
  # ### Starting values and priors
  # ###
  
  # a0 <- start$a0
  # z <- start$z
  # sigma <- start$sigma
  # sigma.mu <- start$sigma.mu
  # pie <- start$pie
  # H <- priors$H
  # mu <- start$mu
  # p <- start$p

 # #Starting values for p and z
  # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  # # p <- sum(idx)/T
  # # z <- numeric(T)
  # # z[idx] <- 1

  
  # ###
  # ###  Setup Variables 
  # ###
  
  # # browser() 
  # T <- nrow(s)  # number of observations
  # mu.0 <- unique(h)  # unique cluster locations
# # browser()
  # n.cls <- nrow(mu.0)  # number of clusters
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  # tab.cls <- table(h.idx)  # tabulate cluster membership
 
  # # Order by decreasing membership
  # ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
  # tab.cls <- tab.cls[ord]
  # idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
 
  # # Propose values for mu.0
  # n.cls.star <- H-n.cls  # number of new clusters to propose
  # mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star
  # samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership  

  # ###
  # ### Create receptacles for output
  # ###
  
  # h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
# h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  # mu.save <- array(0,dim=c(T,2,n.mcmc))  # true animal locations
  # a0.save <- numeric(n.mcmc)  # concentration parameter
  # sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  # sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  # p.save <- numeric(n.mcmc)  # probability of being hauled-out
  # z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  # mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  
  # keep <- list(mu=0,sigma=0,sigma.mu=0)

    
  # ###
  # ### Begin MCMC loop
  # ###
  
  # for (k in 1:n.mcmc) {
    # if(k%%1000==0) cat(k,"");flush.console()

	# ###
	# ### Dirichlet process parameters
	# ###

	# # Update follows the blocked Gibbs sampler of Ishwaran and James (2001) and
    # # Gelman et al. (2014), Section 23.3

	# # Note: sampling order matters here. Clusters must be sampled in same 
	# # order as pie, i.e., sorted by decreasing membership

	   	# ### Sample h.t (cluster assignment indicator)
# # browser()
	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # # (dnorm(mu[x,1],mu.0[samp.cls,1],sigma.mu)*
		  # # dnorm(mu[x,2],mu.0[samp.cls,2],sigma.mu))^(1-z[x])))

		# # Sampled with truncated normal density
	
	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]))

	   	# h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # (dtnorm(rep(mu[x,1],H),mu.0[samp.cls,1],sigma.mu,
		  # lower=min(S[,1]),upper=max(S[,1]))*
		  # dtnorm(rep(mu[x,2],H),mu.0[samp.cls,2],sigma.mu,lower=min(S[,2]),
		  # upper=max(S[,2])))^(1-z[x])))
		# tab.cls <- table(h.idx)  # tabulate cluster membership
		# n.cls <- length(tab.cls)  # number of clusters
	
	  	# # Sort in order of decreasing cluster membership
		# ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
		# tab.cls <- tab.cls[ord]
		# idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
		# samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership
  
 	    # ### Sample pie (stick-breaking process)
    
	  	# tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))  # membership in decreasing order
	    # v <- c(rbeta(H-1,1+tab.cls.tmp,a0+T-cumsum(tab.cls.tmp)),1)  # stick-breaking weights
	    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	    # ### Sample a0 (concentration parameter); See Gelman section 23.3
       
    	# a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
# # a0 <- start$a0

    # ###
    # ### Sample mu.0 (true location of occupied clusters)
    # ###
	   
	# # Sampling order does not matter here

	# get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# idx.0 <- which(h.idx==x&z==0)
		# idx.1 <- which(h.idx==x&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# Sigma.inv <- solve(sigma^2*diag(2))
		# Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		# b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(mu[idx.0,]%*%Sigma.mu.inv)
		# A <- n.1*Sigma.inv+n.0*Sigma.mu.inv
		# A.inv <- solve(A)
		# mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# mu.0.tmp
	# }

# # browser()	
	# mu.0.tmp <- t(sapply(idx.cls,function(x) 
		# get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  # proposals for mu.0	

	# # idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		# # mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	

	# # mu.0[idx,] <- mu.0.tmp[idx,]

	# mu.0[idx.cls,] <- mu.0.tmp
	
	# # mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # update mu.0
	# # mu.0[idx.cls,] <- mu.0.tmp  # update mu.0

	# n.cls.star <- H-n.cls  # number of new clusters to propose
	# mu.0[-idx.cls,] <- cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2]))  # update mu.0 with mu.star


    # ###
    # ### Sample mu (true location of individual)
    # ###
    
# # browser()
	# # Sampling order does not matter here
	
    # # Update mu[t] for z[t]==1, locations at haul-out sites

	# idx <- which(z==1)  # locations at haul-out
	# h.idx.tmp <- h.idx[idx]
	# mu[idx,] <- mu.0[h.idx.tmp,]  # mu is haul-out site for hauled-out individuals

    # # Update mu[t] for z[t]==0, at-sea locations
	# # idx <- which(z==0)  # at-sea locations
	# # h.idx.tmp <- h.idx[idx]
    # # b <- s[idx,]%*%solve(sigma^2*diag(2))+mu.0[h.idx.tmp,]%*%solve(sigma.mu^2*diag(2))
    # # A.inv <- solve(solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2)))  # var-cov matrix    
    # # mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))  # mean matrix
	# # T.0 <- nrow(mu.tmp)   
    # # mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	# # rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))  # proposals for mu
    # # idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      # # mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	# # mu[idx[idx.tmp],] <- mu.star[idx.tmp,]
# # mu <- start$mu
    
    # ###
    # ### Sample sigma (observation error)
    # ###
# # browser()
    # sigma.star <- rnorm(1,sigma,tune$sigma)
    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # mh.star.sigma <- sum(dnorm(s[,1],mu[,1],sigma.star,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma.star,log=TRUE))
      # mh.0.sigma <- sum(dnorm(s[,1],mu[,1],sigma,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma,log=TRUE))
      # # mh.star.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma.star^2*diag(2),log=TRUE)))
      # # mh.0.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma^2*diag(2),log=TRUE)))
      # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # sigma <- sigma.star
        # keep$sigma <- keep$sigma+1
      # } 
    # }
# # sigma <- start$sigma
    
    # ###
    # ### Sample sigma.mu (disperson around homerange center)
    # ###
# # browser()
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,log=TRUE))
      # # mh.0.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

    # # Sample with truncated normal density
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

# # sigma.mu <- start$sigma.mu

	# ###
    # ### Sample z (haul-out indicator variable)
    # ###

# # browser()
    # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	# # mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    # # n.tmp <- sum(idx)
    # # z[!idx] <- 0
	# # p.tmp1 <- (p*dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		# # dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE))^z[idx]
	# # # p.tmp2 <- ((1-p)*dnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu)*
		# # # dnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu))^(1-z[idx])
	# # # Sample with truncated normal density
	# # p.tmp2 <- ((1-p)*dtnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu,
		# # lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
		# # dtnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu,
		# # lower=min(S[,2]),upper=max(S[,2]),log=FALSE))^(1-z[idx])
	# # p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
	# # z[idx] <- rbinom(n.tmp,1,p.tmp)

    # # # fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
	    # # # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    # # # p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    # # # z[idx] <- rbinom(n.tmp,1,p.tmp)

# # z <- start$z



	# # ###
    # # ### Sample p (probability of hauled out)
    # # ###
    
    # # sumz <- sum(z)
    # # p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)
# # p <- start$p

    # ###
    # ###  Save samples 
    # ###

	# h.idx.save[,k] <- h.idx 
# h.save[,,k] <- mu.0[h.idx,]
	# mu.0.save[,,k] <- mu.0
    # a0.save[k] <- a0    
    # sigma.save[k] <- sigma
    # sigma.mu.save[k] <- sigma.mu
    # mu.save[,,k] <- mu
    # p.save[k] <- p
	# z.save[,k] <- z
  # }
  
  # ###
  # ### Write output
  # ###
  
  # keep$sigma <- keep$sigma/n.mcmc
  # keep$sigma.mu <- keep$sigma.mu/n.mcmc
  # cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  # cat(paste("\nsigma.mu acceptance rate:",keep$sigma.mu)) 
  # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
  # list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,mu=mu.save,a0=a0.save,sigma=sigma.save,
  	# sigma.mu=sigma.mu.save,p=p.save,z=z.save,keep=keep,n.mcmc=n.mcmc)
  
# }




# haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  # ###
  # ### Libraries and Subroutines
  # ###
  
  # library(MCMCpack)  # for Dirichlet distribution functions
  # library(data.table)  # for tabulating and summing
  # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # library(doParallel)  # for parallel processing
  # library(foreach)  # for parallel processing  
  # # library(mvtnorm)  # for multivariate normal density
  # library(msm)  # for truncated normal density
  
	# get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# idx.0 <- which(h.idx==x&z==0)
		# idx.1 <- which(h.idx==x&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# b <- colSums(s[idx.1,]%*%solve(sigma^2*diag(2)))+
			# colSums(mu[idx.0,]%*%solve(sigma.mu^2*diag(2)))
		# A <- n.1*solve(sigma^2*diag(2))+n.0*solve(sigma.mu^2*diag(2))
		# A.inv <- solve(A)
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# mu.0.tmp
	# }
  
  
  # ###
  # ###  Create cluster for parallel processing
  # ###
  
  # #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  # #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  # #   mcoptions <- list(preschedule=TRUE)
  # #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # ###
  # ### Starting values and priors
  # ###
  
  # a0 <- start$a0
  # # h <- start$h
  # z <- start$z
  # sigma <- start$sigma
  # sigma.mu <- start$sigma.mu
  # pie <- start$pie
  # H <- priors$H
  # mu <- start$mu
  # p <- start$p

 # #Starting values for p and z
  # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  # # p <- sum(idx)/T
  # # z <- numeric(T)
  # # z[idx] <- 1

  
  # ###
  # ###  Setup Variables 
  # ###
  
  # # browser() 
  # T <- nrow(s)  # number of observations
  # mu.0 <- unique(h)  # unique cluster locations
# # browser()
  # n.cls <- nrow(mu.0)  # number of clusters
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  # tab.cls <- table(h.idx)  # tabulate cluster membership
 
  # # Order by decreasing membership
  # ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
  # mu.0 <- mu.0[ord,]
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  # tab.cls <- table(h.idx)  # tabulate cluster membership
  
  # # Add proposals for mu.0	
  # n.cls.star <- H-n.cls  # number of new clusters to propose
  # mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star

  # # idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
  # # samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership
  
  
  # ###
  # ### Create receptacles for output
  # ###
  
  # h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
# h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  # mu.save <- array(0,dim=c(T,2,n.mcmc))  # true animal locations
  # a0.save <- numeric(n.mcmc)  # concentration parameter
  # sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  # sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  # p.save <- numeric(n.mcmc)  # probability of being hauled-out
  # z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  # mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  
  # keep <- list(mu=0,sigma=0,sigma.mu=0)

    
  # ###
  # ### Begin MCMC loop
  # ###
  
  # for (k in 1:n.mcmc) {
    # if(k%%1000==0) cat(k,"");flush.console()

	# ###
	# ### Dirichlet process parameters
	# ###

	# # Update follows the blocked Gibbs sampler of Ishwaran and James (2001) and
    # # Gelman et al. (2014), Section 23.3

	# # Note: sampling order matters here. Clusters must be sampled in same 
	# # order as pie, i.e., sorted by decreasing membership

	   	# ### Sample h.t (cluster assignment indicator)
# # browser()
	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # # (dnorm(mu[x,1],mu.0[samp.cls,1],sigma.mu)*
		  # # dnorm(mu[x,2],mu.0[samp.cls,2],sigma.mu))^(1-z[x])))

		# # Sampled with truncated normal density
	   	# h.idx <- sapply(1:T,function(x) sample(1:H,1,
	      # prob=pie*(dnorm(s[x,1],mu.0[,1],sigma)*
	      # dnorm(s[x,2],mu.0[,2],sigma))^z[x]))

	   	# # h.idx <- sapply(1:T,function(x) sample(samp.cls,1,
	      # # prob=pie*(dnorm(s[x,1],mu.0[samp.cls,1],sigma)*
	      # # dnorm(s[x,2],mu.0[samp.cls,2],sigma))^z[x]*
		  # # (dtnorm(rep(mu[x,1],H),mu.0[samp.cls,1],sigma.mu,
		  # # lower=min(S[,1]),upper=max(S[,1]))*
		  # # dtnorm(rep(mu[x,2],H),mu.0[samp.cls,2],sigma.mu,lower=min(S[,2]),
		  # # upper=max(S[,2])))^(1-z[x])))
		# tab.cls <- table(h.idx)  # tabulate cluster membership
		# n.cls <- length(tab.cls)  # number of clusters
	
	  	# # Sort in order of decreasing cluster membership
		# ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
		# tab.cls <- tab.cls[ord]
		# idx.cls <- as.numeric(names(tab.cls))  # 'occupied' clusters
		# mu.0 <- mu.0[idx.cls,]
		# h.idx <- c(1:n.cls)[match(h.idx,idx.cls)]  # cluster membership indicator

		# # samp.cls <- c(idx.cls,setdiff(1:H,idx.cls))  # order of decreasing membership
  
 	    # ### Sample pie (stick-breaking process)
    
	  	# tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))  # membership in decreasing order
	    # v <- c(rbeta(H-1,1+tab.cls.tmp,a0+T-cumsum(tab.cls.tmp)),1)  # stick-breaking weights
	    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	    # ### Sample a0 (concentration parameter); See Gelman section 23.3
       
    	# a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
# # a0 <- start$a0

    # ###
    # ### Sample mu.0 (true location of occupied clusters)
    # ###
	   
	# # Sampling order does not matter here

	# get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		# idx.0 <- which(h.idx==x&z==0)
		# idx.1 <- which(h.idx==x&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# Sigma.inv <- solve(sigma^2*diag(2))
		# Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		# b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(mu[idx.0,]%*%Sigma.mu.inv)
		# A <- n.1*Sigma.inv+n.0*Sigma.mu.inv
		# A.inv <- solve(A)
		# mu.0.tmp <- rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
		# # mu.0.tmp <- rnorm(2,A.inv%*%b,sigma/sqrt(n.1))	# proposal for mu.0	
		# mu.0.tmp
	# }

# # browser()	
	# mu.0.tmp <- t(sapply(1:n.cls,function(x) 
		# get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  # proposals for mu.0	

	# # mu.0.tmp <- t(sapply(idx.cls,function(x) 
		# # get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  # proposals for mu.0	
# # points(mu.0.tmp,pch=as.character(1:4),col=2)
	# idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		# mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	
	# mu.0[idx,] <- mu.0.tmp[idx,]
	
	# # mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # update mu.0
	# # mu.0[idx.cls,] <- mu.0.tmp  # update mu.0

	# n.cls.star <- H-n.cls  # number of new clusters to propose
	# mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star

	# # mu.0 <- rbind(mu.0[idx.cls,], cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star

    # ###
    # ### Sample mu (true location of individual)
    # ###
    
# # browser()
	# # Sampling order does not matter here
	
    # # Update mu[t] for z[t]==1, locations at haul-out sites
	# idx <- which(z==1)  # locations at haul-out
	# h.idx.tmp <- h.idx[idx]
	# mu[idx,] <- mu.0[h.idx.tmp,]  # mu is haul-out site for hauled-out individuals

    # # Update mu[t] for z[t]==0, at-sea locations
	# # idx <- which(z==0)  # at-sea locations
	# # h.idx.tmp <- h.idx[idx]
    # # b <- s[idx,]%*%solve(sigma^2*diag(2))+mu.0[h.idx.tmp,]%*%solve(sigma.mu^2*diag(2))
    # # A.inv <- solve(solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2)))  # var-cov matrix    
    # # mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))  # mean matrix
	# # T.0 <- nrow(mu.tmp)   
    # # mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	# # rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))  # proposals for mu
    # # idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      # # mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	# # mu[idx[idx.tmp],] <- mu.star[idx.tmp,]
# # mu <- start$mu
    
    # ###
    # ### Sample sigma (observation error)
    # ###
# # browser()
    # sigma.star <- rnorm(1,sigma,tune$sigma)
    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # mh.star.sigma <- sum(dnorm(s[,1],mu[,1],sigma.star,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma.star,log=TRUE))
      # mh.0.sigma <- sum(dnorm(s[,1],mu[,1],sigma,log=TRUE)+
	      # dnorm(s[,2],mu[,2],sigma,log=TRUE))
      # # mh.star.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma.star^2*diag(2),log=TRUE)))
      # # mh.0.sigma <- sum(sapply(1:T, function(x) 
      	# # dmvnorm(s[x,],mu[x,],sigma^2*diag(2),log=TRUE)))
      # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # sigma <- sigma.star
        # keep$sigma <- keep$sigma+1
      # } 
    # }
# # sigma <- start$sigma
    
    # ###
    # ### Sample sigma.mu (disperson around homerange center)
    # ###
# # browser()
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,log=TRUE))
      # # mh.0.sigma.mu <- sum(dnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,log=TRUE)+
	      # # dnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

    # # Sample with truncated normal density
    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
	  # # idx <- which(z==0)
	  # # mu.0.tmp <- mu.0[h.idx[idx],]
      # # mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu.star,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu.star,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0.tmp[,1],sigma.mu,
	      # # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	      # # dtnorm(mu[idx,2],mu.0.tmp[,2],sigma.mu,
	      # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
      # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        # # sigma.mu <- sigma.mu.star
        # # keep$sigma.mu <- keep$sigma.mu+1
      # # } 
    # # }

# # sigma.mu <- start$sigma.mu

	# ###
    # ### Sample z (haul-out indicator variable)
    # ###

# # browser()
    # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	# # mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    # # n.tmp <- sum(idx)
    # # z[!idx] <- 0
	# # p.tmp1 <- (p*dnorm(s[idx,1],mu[idx,1],sigma,log=FALSE)*
		# # dnorm(s[idx,2],mu[idx,2],sigma,log=FALSE))^z[idx]
	# # # p.tmp2 <- ((1-p)*dnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu)*
		# # # dnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu))^(1-z[idx])
	# # # Sample with truncated normal density
	# # p.tmp2 <- ((1-p)*dtnorm(mu[idx,1],mu.0[h.idx[idx],1],sigma.mu,
		# # lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
		# # dtnorm(mu[idx,2],mu.0[h.idx[idx],2],sigma.mu,
		# # lower=min(S[,2]),upper=max(S[,2]),log=FALSE))^(1-z[idx])
	# # p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
	# # z[idx] <- rbinom(n.tmp,1,p.tmp)

    # # # fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
	    # # # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    # # # p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    # # # z[idx] <- rbinom(n.tmp,1,p.tmp)

# # z <- start$z



	# # ###
    # # ### Sample p (probability of hauled out)
    # # ###
    
    # # sumz <- sum(z)
    # # p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)
# # p <- start$p

    # ###
    # ###  Save samples 
    # ###

	# h.idx.save[,k] <- h.idx 
# h.save[,,k] <- mu.0[h.idx,]
	# mu.0.save[,,k] <- mu.0
    # a0.save[k] <- a0    
    # sigma.save[k] <- sigma
    # sigma.mu.save[k] <- sigma.mu
    # mu.save[,,k] <- mu
    # p.save[k] <- p
	# z.save[,k] <- z
  # }
  
  # ###
  # ### Write output
  # ###
  
  # keep$sigma <- keep$sigma/n.mcmc
  # keep$sigma.mu <- keep$sigma.mu/n.mcmc
  # cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  # cat(paste("\nsigma.mu acceptance rate:",keep$sigma.mu)) 
  # list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,mu=mu.save,a0=a0.save,sigma=sigma.save,
  	# sigma.mu=sigma.mu.save,p=p.save,z=z.save,keep=keep,n.mcmc=n.mcmc)
  
# }












# haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  # ###
  # ### Libraries and Subroutines
  # ###
  
  # library(MCMCpack)  # for Dirichlet distribution functions
  # library(data.table)  # for tabulating and summing
  # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  # library(doParallel) #For parallel processing
  # library(foreach) #For parallel processing  
  
  # # get.mu.0 <- function(x,dt,tab,sigma){
    # # y.tmp <- dt[z1==tab[x,z1]&z2==tab[x,z2],.(y1,y2)]
    # # y.tmp <- as.matrix(y.tmp)
    # # A <- solve(sigma^2*diag(2))*tab[x,N]
    # # b <- colSums(y.tmp%*%solve(sigma^2*diag(2)))
    # # t(solve(A)%*%b)
  # # }
  
  
  # ###
  # ###  Create cluster for parallel processing
  # ###
  
  # #   if(is.null(n.cores)) n.cores <- detectCores() - 1
  # #   if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore functionality	
  # #   mcoptions <- list(preschedule=TRUE)
  # #   cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
  
  # ###
  # ### Starting values and priors
  # ###
  
  # a0 <- start$a0
  # # h <- start$h
  # z <- start$z
  # sigma <- start$sigma
  # sigma.mu <- start$sigma.mu
  # pie <- start$pie
  # H <- priors$H
  # mu <- start$mu
  # p <- start$p

 # #Starting values for p and z
  # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
  # # p <- sum(idx)/T
  # # z <- numeric(T)
  # # z[idx] <- 1

  
  # ###
  # ###  Setup Variables 
  # ###
  
  # # browser() 
  # T <- nrow(s)  # number of observations
  # mu.0 <- unique(h)  # unique cluster locations
  # n.cls <- nrow(mu.0)	
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]
  # tab.cls <- table(h.idx)	
  # ord <- order(tab.cls,decreasing=TRUE)	

  # # dt <- as.data.table(cbind(sx=s[,1],sy=s[,2],mx=mu[,1],my=mu[,2],hx=h[,1],hy=h[,2],z=z))
  # # tab.z <- dt[,.N,by=.(hx,hy,z)]  
  # # tab <- tab.z[,.(N=sum(N)),by=.(hx,hy)]  
  # # setkey(tab,N)
  # # n.cluster <- tab[,.N]  
  # # ord <- n.cluster:1
  
  # h.idx.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
  # mu.save <- array(0,dim=c(T,2,n.mcmc))  # true animal locations
  # a0.save <- numeric(n.mcmc)  # concentration parameter
  # sigma.save <- numeric(n.mcmc)  # telemetry measurement error
  # sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
  # p.save <- numeric(n.mcmc)  # probability of being hauled-out
  # z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
  # mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
  
  # keep <- list(mu=0,sigma=0,sigma.mu=0)

    
  # ###
  # ### Begin MCMC loop
  # ###
  
  # for (k in 1:n.mcmc) {
    # if(k%%1000==0) cat(k,"");flush.console()
    
    # ###
    # ### Sample mu (true location)
    # ###
    
    # # browser()
    # # Update mu[t] for z[t]==1, locations at haul-out sites
	# idx <- which(z==1)  # locations at haul-out
	# mu[idx,] <- h[idx,]

	# # dt[idx,mx:=hx]  # just the haul-out location
 	# # dt[idx,my:=hy]  # just the haul-out location

    # # Update mu[t] for z[t]==0, locations at-sea
	# idx <- which(z==0)
    # b <- s[idx,]%*%solve(sigma^2*diag(2))+h[idx,]%*%solve(sigma.mu^2*diag(2))
    # # b <- s[idx,]%*%solve(sigma^2*diag(2))+
    	# # as.matrix(dt[idx,.(hx,hy)],T,2,byrow=TRUE)%*%solve(sigma.mu^2*diag(2))
    # A <- solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2))    
    # A.inv <- solve(A)
        
    # mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))
	# T.0 <- nrow(mu.tmp)   
    # mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	# rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))
    # idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      # mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	# mu[idx[idx.tmp],] <- mu.star[idx.tmp,]

	# # dt[idx,mx:=mu.star[idx,1]]
 	# # dt[idx,my:=mu.star[idx,2]]
 
    # # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
    
    # ###
    # ### Sample theta (cluster parameter)
    # ###
    
    # # Sampler currently disregards support P0
# # browser()

  # # get.mu.0 <- function(x,tab,tab.z,dt,sigma,sigma.mu){
	# # browser()
	# # s.tmp <- as.matrix(dt[hx==tab[x,hx]&hy==tab[x,hy]&z==1,.(sx,sy)])
	# # mu.tmp <- as.matrix(dt[hx==tab[x,hx]&hy==tab[x,hy]&z==0,.(mx,my)])
	# # n.0 <- nrow(s.tmp)
	# # n.1 <- nrow(mu.tmp)
	# # # n.0 <- tab.z[hx==tab[x,hx]&z==0,N]
	# # # n.0 <- ifelse(length(n.0)>0,n.0,0)
	# # # n.1 <- tab.z[hx==tab[x,hx]&z==1,N]
	# # # n.1 <- ifelse(length(n.1)>0,n.1,0)
    # # A <- solve(sigma^2*diag(2))*n.1+solve(sigma.mu^2*diag(2))*n.0
	# # b <- colSums(s.tmp%*%solve(sigma^2*diag(2)))+colSums(mu.tmp%*%solve(sigma.mu^2*diag(2)))
	# # A.inv <- solve(A)
	# # mu.tmp <- A.inv%*%b
	# # c(rnorm(1,mu.tmp[1],sqrt(A.inv[1,1])),rnorm(1,mu.tmp[2],sqrt(A.inv[2,2])))
  # # }

	# mu.0.tmp <- unique(h)
	# n.tmp <- numeric(nrow(mu.0.tmp))
	# for(i in 1:nrow(mu.0.tmp)){
		# idx.0 <- which(h[,1]==mu.0.tmp[i,1]&h[,2]==mu.0.tmp[i,2]&z==0)
		# idx.1 <- which(h[,1]==mu.0.tmp[i,1]&h[,2]==mu.0.tmp[i,2]&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# b <- colSums(s[idx.1,]%*%solve(sigma^2*diag(2)))+
			# colSums(mu[idx.0,]%*%solve(sigma.mu^2*diag(2)))
		# A.inv <- solve(n.1*solve(sigma^2*diag(2))+n.0*solve(sigma.mu^2*diag(2)))
		# tmp <- rnorm(2,A.inv%*%b,diag(A.inv))		
		# if(tmp[1]>S.tilde[1,1]&tmp[1]<S.tilde[2,1]&
			# tmp[2]>S.tilde[1,2]&tmp[2]<S.tilde[3,2]){
				# mu.0.tmp[i,] <- tmp				
		# }
		# n.tmp[i] <- n.0+n.1
	# }
	
	# ord <- order(n.tmp,decreasing=TRUE)
    # mu.0.tmp <- mu.0.tmp[ord,] 
    # n.tmp <- n.tmp[ord]
    # n.cls <- length(n.tmp)
    
    # # mu.0 <- t(sapply(1:n.cluster,function(x) get.mu.0(x,tab,tab.z,dt,sigma,sigma.mu)))
 	# # mu.0 <- mu.0[ord,]
    # n.cls.star <- H-n.cls
    # mu.0 <- rbind(mu.0.tmp,cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))
  
  

    # ###
    # ### Sample h.t (cluster assignments)
    # ###
     
    # # browser()
   	# h.t <- sapply(1:T,function(x) sample(1:H,1,
      # prob=pie*(dnorm(s[x,1],mu.0[,1],sigma)*dnorm(s[x,2],mu.0[,2],sigma))^z[x]*
	  # (dnorm(mu[x,1],mu.0[,1],sigma.mu)*dnorm(mu[x,2],mu.0[,2],sigma.mu))^(1-z[x])))
	# h <- mu.0[h.t,]

   	# # h.t <- sapply(1:T,function(x) sample(1:H,1,
      # # prob=pie*(dnorm(s[x,1],mu.0[,1],sigma)*dnorm(s[x,2],mu.0[,2],sigma))^z[x]*
      # # (dnorm(dt[x,mx],mu.0[,1],sigma.mu)*dnorm(dt[x,my],mu.0[,2],sigma.mu))^z[x]))
    # # dt[,hx:=mu.0[h.t,1]]
    # # dt[,hy:=mu.0[h.t,2]]    


    # ###
    # ### Sample pie (stick-breaking process)
    # ###

   # # tab.z <- dt[,.N,by=.(hx,hy,z)]  
   # # tab <- tab.z[,.(N=sum(N)),by=.(hx,hy)]  
   # # setkey(tab,N)
   # # n.cluster <- tab[,.N]  
   # # ord <- n.cluster:1
   # # tab.tmp <- c(tab[ord,N],rep(0,H-n.cluster-1))

    # # Update stick-breaking weights
    
    # # v <- c(rbeta(H-1,1+tab.tmp,a0+T-cumsum(tab.tmp)),1)
    # # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	# mu.0.tmp <- unique(h)
	# n.cls <- nrow(mu.0.tmp)
	# n.tmp <- table(h.t)
	# n.tmp <- sort(n.tmp,decreasing=TRUE)
    # n.tmp <- c(n.tmp,rep(0,H-n.cls-1))
    # v <- c(rbeta(H-1,1+n.tmp,a0+T-cumsum(n.tmp)),1)
    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities



    # ###
    # ### Sample a0 (concentration parameter); See Gelman section 23.3
    # ###
    
    # a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
	# # a0 <- start$a0


    # ###
    # ### Sample sigma (observation error)
    # ###

    # # sigma.star <- rnorm(1,sigma,tune$sigma)
    # # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # # mh.star.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma.star,log=TRUE)+
        # # dnorm(dt[,y2],dt[,z2],sigma.star,log=TRUE))
      # # mh.0.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma,log=TRUE)+
        # # dnorm(dt[,y2],dt[,z2],sigma,log=TRUE))
      # # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # # sigma <- sigma.star
        # # keep$sigma <- keep$sigma+1
      # # } 
    # # }
    
    
    # ###
    # ### Sample z (haul-out indicator variable)
    # ###

# #     browser()
    # # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	# # mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    # # n.tmp <- sum(idx)
    # # z[!idx] <- 0
    # # fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
	    # # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    # # p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    # # z[idx] <- rbinom(n.tmp,1,p.tmp)


    # ###
    # ### Sample sigma.mu (disperson around homerange center)
    # ###

    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # # if(sigma.mu.star>0){
		# # mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu.star,
			# # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	        # # dtnorm(mu[idx,2],mu.0[2],sigma.mu.star,
	        # # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
    	# # mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu,
      		# # lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
      		# # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
	    # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
    	    # # sigma.mu <- sigma.mu.star
        	# # keep$sigma.mu <- keep$sigma.mu+1
      	# # } 
	# # }    


	# ###
    # ### Sample p (probability of hauled out)
    # ###
    
    # # sumz <- sum(z)
    # # p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)


    # ###
    # ###  Save samples 
    # ###

# h.save[,,k] <- h 
    # # h.save[,1,k] <- dt[,hx]
    # # h.save[,2,k] <- dt[,hy]
    # a0.save[k] <- a0    
    # sigma.save[k] <- sigma
    # sigma.mu.save[k] <- sigma.mu
    # mu.save[,,k] <- mu
    # p.save[k] <- p
	# z.save[,k] <- z
  # }
  
  # ###
  # ### Write output
  # ###
  
  # keep$sigma <- keep$sigma/n.mcmc
  # cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  # list(h=h.save,mu=mu.save,a0=a0.save,sigma=sigma.save,sigma.mu=sigma.mu.save,p=p.save,
	# z=z.save,keep=keep,n.mcmc=n.mcmc)
  
# }