haulout.dpmixture.3.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S,S.tilde,U,sigma.alpha,
	priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (13 AUG 2015)
 	### See haulout.dpmixture.2.sim.R to simulate data according to this model specification,
 	### and haulout.dp.mixture.2.pdf for the model description, model statement, and
 	### full conditional distributions
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s; X.tilde=design matrix containing covariates
 	### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	### W.tilde=basis expansion for y; S=support of true locations/movement process (mu_t);
 	### S.tilde=support Dirichlet process mixture (i.e., haul-out sites); sigma.alpha=
 	### standard deviation of parameter model for 'random' effects
 	###
 
 
	t.start <- Sys.time()
	

	###
	### Libraries and Subroutines
	###
  
	# library(MCMCpack)  # for Dirichlet distribution functions
	library(data.table)  # for tabulating and summing
	# library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	library(doParallel)  # for parallel processing
	library(foreach)  # for parallel processing  
	# library(mvtnorm)  # for multivariate normal density
	# library(msm)  # for truncated normal density
  
	# get.mu.0 <- function(x,h.idx,z,s,Sigma.inv,Sigma.mu.inv){
		# # browser()
		# idx.0 <- which(h.idx==x&z==0)
		# idx.1 <- which(h.idx==x&z==1)
		# n.0 <- length(idx.0)
		# n.1 <- length(idx.1)
		# b <- colSums(s[idx.1,]%*%Sigma.inv)+colSums(s[idx.0,]%*%(Sigma.inv+Sigma.mu.inv))
		# A.inv <- solve(n.1*Sigma.inv+n.0*(Sigma.inv+Sigma.mu.inv))
		# rnorm(2,A.inv%*%b,sqrt(diag(A.inv)))	# proposal for mu.0	
	# }

	truncnormsamp <- function(mu,sig2,low,high,nsamp){
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}
	
	###
	###  Create cluster for parallel processing
	###
  
	# if(is.null(n.cores)) n.cores <- detectCores() - 1
	# if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore 	functionality	
	# mcoptions <- list(preschedule=TRUE)
	# cat(paste("\nUsing",n.cores,"cores for parallel processing."))
  
	###
	###  Setup Variables 
	###
  
	# browser() 
	T <- nrow(s)  # number of telemetry locations
	n <- length(y)  # number of wet/dry observations
	qX <- ncol(X)
	qW <- ncol(W)
	qU <- ncol(U)
	y1 <- (y==1)
	y0 <- (y==0)
	y1.sum <- sum(y1)
	y0.sum <- sum(y0)
	v <- numeric(n)  # auxilliary variable for y
# browser()	
	tune$mu.0 <- seq(-tune$mu.0,tune$mu)
	idx <- which(tune$mu.0==0)
	tune$mu.0 <- tune$mu.0[-idx]

	###
	### Starting values and priors
	###
# browser() 
	beta <- matrix(start$beta,qX)
	alpha <- matrix(0,qW)
	gamma <- matrix(start$gamma,qU)
	# alpha <- start$alpha
	theta <- start$theta
	sigma <- start$sigma
	sigma.mu <- start$sigma.mu
	pie <- start$pie

	H <- priors$H
	Sigma.inv <- solve(sigma^2*diag(2))
	Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
  	mu.beta <- matrix(0,qX,1)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	mu.alpha <- matrix(0,qW,1)
	Sigma.alpha <- diag(qW)*sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)
	mu.gamma <- matrix(priors$mu.gamma,qU,1)
	Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	Sigma.gamma.inv <- solve(Sigma.gamma)


	###
	### Set up Dirichlet process mixture variables
	###
# browser()
	mu.0 <- unique(start$h)  # unique cluster location s
	n.cls <- nrow(mu.0)  # number of clusters
	h.idx <- c(1:n.cls)[match(start$h[,2],mu.0[,2])]  # cluster membership indicator
	S.tilde.idx <- 1:nrow(S.tilde)
	cell.idx <- sapply(mu.0[,2],function(x) which(S.tilde[,2]==x))

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
	idx <- sample(S.tilde.idx[-cell.idx],n.cls.star,replace=FALSE)
	mu.0 <- rbind(mu.0, S.tilde[idx,])
# sum(duplicated(mu.0))
	# mu.0 <- rbind(mu.0, cbind(runif(n.cls.star,min(S.tilde[,1]),max(S.tilde[,1])),
      # runif(n.cls.star,min(S.tilde[,2]),max(S.tilde[,2]))))  # update mu.0 with mu.star

  	###
	### Create receptacles for output
	###
  
	h.idx.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
	h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
	theta.save <- numeric(n.mcmc)  # concentration parameter
	mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
	n.cls.save <- numeric(n.mcmc)
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	alpha.save <- matrix(0,n.mcmc,qW)
	beta.save <- matrix(0,n.mcmc,qX)
	gamma.save <- matrix(0,n.mcmc,qU)
	v.save <- matrix(0,n,n.mcmc)
	sigma.save <- numeric(n.mcmc)  # telemetry measurement error
	sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	keep <- list(sigma=0,sigma.mu=0,mu.0=0)
    
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) cat(k,"");flush.console()

		###
		### Updates pertaining to wet/dry status of y and z
		###

			###  Sample v(t.tilde) (auxilliary variable for y) 
	
			linpred <- X.tilde%*%beta+W.tilde%*%alpha
		  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
		  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	
			###  Sample alpha (coefficients on basis expansion W.tilde)
	
			A.inv <- solve(t(W.tilde)%*%W.tilde+Sigma.alpha.inv)
			b <- t(W.tilde)%*%(v-X.tilde%*%beta)
			alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
		
			###  Sample beta (coefficients on covariates influencing P(y==1))
	
		 	A.inv <- solve(t(X.tilde)%*%X.tilde+Sigma.beta.inv)
		  	b <- t(X.tilde)%*%(v-W.tilde%*%alpha)  # +mu.beta%*%Sigma.beta.inv
		  	beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
	
		    ### Sample z(t) (prediction of haul-out indicator variable for times t)
	
			p <- pnorm(X%*%beta+W%*%alpha)
			p.tmp1 <- p*dnorm(s[,1],mu.0[h.idx,1],sigma,log=FALSE)*
				dnorm(s[,2],mu.0[h.idx,2],sigma,log=FALSE)
			p.tmp2 <- (1-p)*dnorm(s[,1],mu.0[h.idx,1],sqrt(sigma^2+sigma.mu^2),log=FALSE)*
				dnorm(s[,2],mu.0[h.idx,2],sqrt(sigma^2+sigma.mu^2),log=FALSE)
			p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
			z <- rbinom(T,1,p.tmp)
# z <- start$z
		
		###
		### Dirichlet process parameters
		### Note: sampling order matters here. Cluster parameters must be 
		### sampled in same order as pie, i.e., sorted by decreasing membership
		###
	
			# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
			# and Gelman et al. (2014), Section 23.3
			
		   	### Sample h.t (cluster assignment indicator)
# browser()
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
			
		    eta <- c(rbeta(H-1,1+tab.cls.tmp,theta+T-cumsum(tab.cls.tmp)),1)  # stick-breaking
		    	# weights
		    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities
	
		    ### Sample theta (concentration parameter); See Gelman section 23.3
	       
	    	theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta

		###
	    ### Sample mu.0 (true location of occupied clusters)
		### Note: sampling order does not matter here
	    ###
		   
		# Use for data.table functionality
		# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0
			# get.mu.0(x,dt.h.idx[1:T,h.idx],z,s,mu,sigma,sigma.mu,S.tilde)))
# browser()			
	get.mh.mu.0 <- function(x,s,h.idx,idx.cls,mu.0,mu.0.star,z,cell.idx,cell.idx.star,
		sigma,sigma.mu,U,gamma){
	    # browser()
		# x
		idx.0 <- which(h.idx==idx.cls[x]&z==0)
		idx.1 <- which(h.idx==idx.cls[x]&z==1)
    	sd.tmp <- sqrt(sigma^2+sigma.mu^2)
		mh.star <- sum(dnorm(s[idx.1,2],mu.0.star[x,2],sigma,log=TRUE),
			dnorm(s[idx.0,2],mu.0.star[x,2],sd.tmp,log=TRUE)) +			
			t(gamma)%*%U[cell.idx.star[x],] - sum(U%*%gamma)
		mh.0 <- sum(dnorm(s[idx.1,2],mu.0[x,2],sigma,log=TRUE),
			dnorm(s[idx.0,2],mu.0[x,2],sd.tmp,log=TRUE)) +			
			t(gamma)%*%U[cell.idx[x],] - sum(U%*%gamma)
		exp(mh.star-mh.0)>runif(1)
	}

		# Use for base functionality	
		# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0	
			# get.mu.0(x,h.idx,z,s,Sigma.inv,Sigma.mu.inv)))  

		cell.idx <- sapply(mu.0[idx.cls,2],function(x) which(S.tilde[,2]==x))
		mu.0.star <- matrix(NA,n.cls,2)  # in same order as idx.cls
		cell.idx.star <- cell.idx+sample(tune$mu.0,n.cls,replace=TRUE)  
		idx <- which(cell.idx.star>0&cell.idx.star<max(S.tilde.idx)&
			!duplicated(cell.idx.star))  # idx of unique mu.0.star in S.tilde	
		mu.0.star[idx,] <- S.tilde[cell.idx.star[idx],]
		mu.0.tmp <- mu.0[idx.cls,]
		mh.mu.0 <- t(sapply(idx,function(x)  # proposals for mu.0	
			get.mh.mu.0(x,s,h.idx,idx.cls,mu.0.tmp,mu.0.star,z,cell.idx,cell.idx.star,
			sigma,sigma.mu,U,gamma)))   
		idx <- idx[mh.mu.0]
		n.keep <- sum(mh.mu.0)
		# if(is.na(n.keep)) browser()
		# print(n.keep)
		# if(n.keep>0){
			keep$mu.0 <- keep$mu.0+n.keep
			mu.0[idx.cls[idx],] <- mu.0.star[idx,]  # accept proposals in S.tilde
			cell.idx[idx] <- cell.idx.star[idx]

			# Propose new values for unoccupied clusters
			n.cls.star <- H-n.cls  # number of new clusters to propose
			idx <- sample(S.tilde.idx[-cell.idx],n.cls.star,replace=FALSE)
			idx <- sample(setdiff(S.tilde.idx,cell.idx),n.cls.star,replace=FALSE)
			mu.0[-idx.cls,] <- S.tilde[idx,]
		# }
print(c(n.keep,sum(duplicated(mu.0))))

		###
	    ### Sample sigma (observation error)
	    ###
# browser()
	    sigma.star <- rnorm(1,sigma,tune$sigma)
	    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
			idx <- z==1
			mu.0.tmp <- mu.0[h.idx,]
	    	sd.tmp.star <- sqrt(sigma.star^2+sigma.mu^2)
	    	sd.tmp <- sqrt(sigma^2+sigma.mu^2)
	    	mh.star.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma.star,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma.star,log=TRUE))#+
		    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sd.tmp.star,log=TRUE)+
		    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sd.tmp.star,log=TRUE))
		    mh.0.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma,log=TRUE))#+
		    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sd.tmp,log=TRUE)+
		    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sd.tmp,log=TRUE))
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
	    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    if(sigma.mu.star>priors$sigma.mu.l & sigma.star<priors$sigma.mu.u){
			idx <- which(z==0)
			mu.0.tmp <- mu.0[h.idx[idx],]
	    	sd.tmp.star <- sqrt(sigma^2+sigma.mu.star^2)
	    	sd.tmp <- sqrt(sigma^2+sigma.mu^2)
		    mh.star.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sd.tmp.star,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sd.tmp.star,log=TRUE))
		    mh.0.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sd.tmp,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sd.tmp,log=TRUE))
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
	        	Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		        keep$sigma.mu <- keep$sigma.mu+1
	    	} 
	    }
# sigma.mu <- start$sigma.mu
	
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
		theta.save[k] <- theta    
		sigma.save[k] <- sigma
		sigma.mu.save[k] <- sigma.mu
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		gamma.save[k,] <- gamma
		v.save[,k] <- v
		z.save[,k] <- z
		n.cls.save[k] <- n.cls
	}
  
	###
	### Write output
	###

	keep$sigma <- keep$sigma/n.mcmc
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu.0 <- keep$mu.0/sum(n.cls.save)
	cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	list(h.idx=h.idx.save,h=h.save,mu.0=mu.0.save,beta=beta.save,alpha=alpha.save,
		gamma=gamma.save,theta=theta.save,sigma=sigma.save,sigma.mu=sigma.mu.save,z=z.save,
		v=v.save,n.cls=n.cls.save,keep=keep,n.mcmc=n.mcmc)
}