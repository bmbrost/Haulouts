haulout.dpmixture.mcmc <- function(s,S.tilde,S,priors,tune,start,n.mcmc,n.cores=NULL){
  
  ###
  ### Libraries and Subroutines
  ###
  
  library(MCMCpack)  # for Dirichlet distribution functions
  library(data.table)  # for tabulating and summing
  library(dplyr)  # dense_rank() for ranking clusters smallest to largest
  library(doParallel) #For parallel processing
  library(foreach) #For parallel processing  
  
  # get.mu.0 <- function(x,dt,tab,sigma){
    # y.tmp <- dt[z1==tab[x,z1]&z2==tab[x,z2],.(y1,y2)]
    # y.tmp <- as.matrix(y.tmp)
    # A <- solve(sigma^2*diag(2))*tab[x,N]
    # b <- colSums(y.tmp%*%solve(sigma^2*diag(2)))
    # t(solve(A)%*%b)
  # }
  
  
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
  # h <- start$h
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
  
  browser() 
  T <- nrow(s)  # number of observations
  mu.0 <- unique(h)  # unique cluster locations
  n.cls <- nrow(mu.0)  # number of clusters
  h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  mu.0 <- rbind(mu.0,matrix(0,H-n.cls,2))  # fill in with empty clusters	
  tab.cls <- table(h.idx)	
  ord <- order(tab.cls,decreasing=TRUE)	 # sort clusters by membership
  tab.cls <- tab.cls[ord]	
  
  # mu.0 <- rbind(mu.0[ord,],matrix(0,H-n.cls,2))  # fill in with empty clusters
  # h.idx <- c(1:n.cls)[match(start$h[,1],mu.0[,1])]  # cluster membership indicator
  	# # 1 is most populated cluster, 2 is second most population cluster, etc.


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
    ### Sample mu (true location)
    ###
    
    # browser()
    # Update mu[t] for z[t]==1, locations at haul-out sites
	idx <- which(z==1)  # locations at haul-out
	mu[idx,] <- mu.0[h.idx[idx],]  # mu is haul-out site for hauled-out individuals

    # Update mu[t] for z[t]==0, at-sea locations
	idx <- which(z==0)  # at-sea locations
    b <- s[idx,]%*%solve(sigma^2*diag(2))+mu.0[h.idx[idx],]%*%solve(sigma.mu^2*diag(2))
    A.inv <- solve(solve(sigma^2*diag(2))+solve(sigma.mu^2*diag(2)))  # var-cov matrix    
    mu.tmp <- t(apply(b,1,function(x) x%*%A.inv))  # mean matrix
	T.0 <- nrow(mu.tmp)   
    mu.star <- cbind(rnorm(T.0,mu.tmp[,1],sqrt(A.inv[1,1])),
    	rnorm(T.0,mu.tmp[,2],sqrt(A.inv[2,2])))  # proposals for mu
    idx.tmp <- which(mu.star[,1]>S[1,1]&mu.star[,1]<S[2,1]&
      mu.star[,2]>S[1,2]&mu.star[,2]<S[3,2])  # mu.star in S
	mu[idx[idx.tmp],] <- mu.star[idx.tmp,]


    # Following Ishwaran and James (2001); Also Gelman et al. (2014), Section 23.3
    
    ###
    ### Sample theta (cluster parameter)
    ###
    
	get.mu.0 <- function(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde){
		idx.0 <- which(h.idx==x&z==0)
		idx.1 <- which(h.idx==x&z==1)
		n.0 <- length(idx.0)
		n.1 <- length(idx.1)
		b <- colSums(s[idx.1,]%*%solve(sigma^2*diag(2)))+
			colSums(mu[idx.0,]%*%solve(sigma.mu^2*diag(2)))
		A <- n.1*solve(sigma^2*diag(2))+n.0*solve(sigma.mu^2*diag(2))
		A.inv <- solve(A)
		mu.0.tmp <- rnorm(2,A.inv%*%b,diag(A.inv))	# proposal for mu.0	
		mu.0.tmp
	}

	
    ###
    ### Order clusters by membership
    ###

	idx.cls <- as.numeric(names(tab.cls))
	
	mu.0.tmp <- t(sapply(idx.cls,function(x) 
		get.mu.0(x,h.idx,z,s,mu,sigma,sigma.mu,S.tilde)))  # proposals for mu.0	
	idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
		mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	
	mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # update mu.0
	
	n.cls.star <- H-n.cls  # number of new clusters to propose
	mu.0 <- rbind(mu.0[idx.cls,], cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      runif(n.cls.star,S.tilde[1,2],S.tilde[3,2])))  # update mu.0 with mu.star

	# mu.0 <- rbind(mu.0[ord,],mu.0[-ord,]) # sort according to membership
      
    ###
    ### Sample h.t (cluster assignments)
    ###
     
    # browser()

   	h.idx <- sapply(1:T,function(x) sample(1:H,1,
      prob=pie*(dnorm(s[x,1],mu.0[,1],sigma)*dnorm(s[x,2],mu.0[,2],sigma))^z[x]*
	  (dnorm(mu[x,1],mu.0[,1],sigma.mu)*dnorm(mu[x,2],mu.0[,2],sigma.mu))^(1-z[x])))
	tab.cls <- table(h.idx)
	n.cls <- length(tab.cls)
	ord <- order(tab.cls,decreasing=TRUE) # sort clusters by membership
	tab.cls <- tab.cls[ord]	
  
    ###
    ### Sample pie (stick-breaking process)
    ###

   # tab.z <- dt[,.N,by=.(hx,hy,z)]  
   # tab <- tab.z[,.(N=sum(N)),by=.(hx,hy)]  
   # setkey(tab,N)
   # n.cluster <- tab[,.N]  
   # ord <- n.cluster:1
   # tab.tmp <- c(tab[ord,N],rep(0,H-n.cluster-1))

    # Update stick-breaking weights
    
    # v <- c(rbeta(H-1,1+tab.tmp,a0+T-cumsum(tab.tmp)),1)
    # pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities

	# mu.0.tmp <- unique(h)
	# n.cls <- nrow(mu.0.tmp)
	# n.tmp <- table(h.t)
	# n.tmp <- sort(n.tmp,decreasing=TRUE)
    # n.tmp <- c(n.tmp,rep(0,H-n.cls-1))

	tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))
    v <- c(rbeta(H-1,1+tab.cls.tmp,a0+T-cumsum(tab.cls.tmp)),1)
    pie <- v*c(1,cumprod((1-v[-H])))  # mixture component probabilities



    ###
    ### Sample a0 (concentration parameter); See Gelman section 23.3
    ###
    
    a0 <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-v[-H])))  
	# a0 <- start$a0


    ###
    ### Sample sigma (observation error)
    ###

    # sigma.star <- rnorm(1,sigma,tune$sigma)
    # if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
      # mh.star.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma.star,log=TRUE)+
        # dnorm(dt[,y2],dt[,z2],sigma.star,log=TRUE))
      # mh.0.sigma <- sum(dnorm(dt[,y1],dt[,z1],sigma,log=TRUE)+
        # dnorm(dt[,y2],dt[,z2],sigma,log=TRUE))
      # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
        # sigma <- sigma.star
        # keep$sigma <- keep$sigma+1
      # } 
    # }
    
    
    ###
    ### Sample z (haul-out indicator variable)
    ###

#     browser()
    # idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
    	# mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
    # n.tmp <- sum(idx)
    # z[!idx] <- 0
    # fmu <- dtnorm(mu[idx,1],mu.0[1],sigma.mu,lower=min(S[,1]),upper=max(S[,1]),log=FALSE)*
	    # dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=FALSE)
    # p.tmp <- (p*fS.tilde)/(p*fS.tilde+(1-p)*fmu)    
    # z[idx] <- rbinom(n.tmp,1,p.tmp)


    ###
    ### Sample sigma.mu (disperson around homerange center)
    ###

    # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
    # if(sigma.mu.star>0){
		# mh.star.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu.star,
			# lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
	        # dtnorm(mu[idx,2],mu.0[2],sigma.mu.star,
	        # lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
    	# mh.0.sigma.mu <- sum(dtnorm(mu[idx,1],mu.0[1],sigma.mu,
      		# lower=min(S[,1]),upper=max(S[,1]),log=TRUE)+
      		# dtnorm(mu[idx,2],mu.0[2],sigma.mu,lower=min(S[,2]),upper=max(S[,2]),log=TRUE))
	    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
    	    # sigma.mu <- sigma.mu.star
        	# keep$sigma.mu <- keep$sigma.mu+1
      	# } 
	# }    


	###
    ### Sample p (probability of hauled out)
    ###
    
    # sumz <- sum(z)
    # p <- rbeta(1,sumz+priors$alpha,T-sumz+priors$beta)


    ###
    ###  Save samples 
    ###

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
  cat(paste("\nsigma acceptance rate:",keep$sigma)) 
  list(h=h.save,mu=mu.save,a0=a0.save,sigma=sigma.save,sigma.mu=sigma.mu.save,p=p.save,
	z=z.save,keep=keep,n.mcmc=n.mcmc)
  
}












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