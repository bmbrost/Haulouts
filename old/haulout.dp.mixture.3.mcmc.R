haulout.dpmixture.3.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S,S.tilde,U,sigma.alpha,
	priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (31 AUG 2015)
 	### See haulout.dpmixture.3.sim.R to simulate data according to this model specification,
 	### and haulout.dp.mixture.3.pdf for the model description, model statement, and
 	### full conditional distributions
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s; X.tilde=design matrix containing covariates
 	### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	### W.tilde=basis expansion for y; S=support of true locations/movement process (mu_t);
 	### S.tilde=support Dirichlet process mixture (i.e., haul-out sites); U=design matrix 
 	###	containing covariates describing all possible haul-out locations; sigma.alpha=
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
  
	get.mh.mu.0 <- function(x,s,z,h.match,cls.idx,mu.0,mu.0.star,sigma,sigma.mu,
		S.match,S.match.star,U,gamma){
		mu.0 <- matrix(mu.0,,2)
		idx.0 <- which(h.match==cls.idx[x]&z==0)
		idx.1 <- which(h.match==cls.idx[x]&z==1)
    	sd.tmp <- sqrt(sigma^2+sigma.mu^2)
		U.0 <- U[S.match[x],]
		U.star <- U[S.match.star[x],]
		mh.star <- sum(dnorm(s[idx.1,2],mu.0.star[x,2],sigma,log=TRUE),
			dnorm(s[idx.0,2],mu.0.star[x,2],sd.tmp,log=TRUE)) +			
			U.star%*%gamma # - int
		mh.0 <- sum(dnorm(s[idx.1,2],mu.0[x,2],sigma,log=TRUE),
			dnorm(s[idx.0,2],mu.0[x,2],sd.tmp,log=TRUE)) +			
			U.0%*%gamma # - int
		exp(mh.star-mh.0)>runif(1)
	}

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
	# if(n.cores==1) registerDoSEQ() else registerDoParallel(cores=n.cores) # multicore 			functionality	
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
	v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	X.comb <- rbind(X,X.tilde)  # combined design matrix for updates on alpha, beta
	W.comb <- rbind(W,W.tilde)  # combined design matrix for updates on alpha, beta 
	S.tilde.idx <- 1:nrow(S.tilde)

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
	z <- start$z
	gamma.int <- log(sum(exp(U%*%gamma)))
	
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
	sd.tmp <- sqrt(sigma^2+sigma.mu^2)

	y1 <- which(y==1)+T
	y0 <- which(y==0)+T
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	linpred <- X.comb%*%beta+W.comb%*%alpha


	###
	### Set up Dirichlet process mixture variables
	###
# browser()
	mu.0 <- unique(start$h)  # clusters
	n.cls <- nrow(mu.0)  # number of clusters
	h.match <- c(1:n.cls)[match(start$h[,2],mu.0[,2])]  # cluster membership indicator
		
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
	cls.tab <- table(h.match)  # tabulate cluster membership
	cls.ord <- order(cls.tab,decreasing=TRUE)  # clusters ordered by membership
	cls.idx <- as.numeric(names(cls.tab))  # idx of occupied clusters
	cls.diff <- setdiff(1:H,cls.idx)  # idx of unoccupied clusters
	cls.samp <- c(cls.idx[cls.ord],cls.diff)  # order in which clusters are sampled
	S.match <- S.tilde.idx[match(mu.0[cls.idx,2],S.tilde[,2])]  # S.tilde membership indicator

 	# Propose values for mu.0
	idx <- sample(S.tilde.idx[-S.match],H-n.cls,replace=FALSE)  # idx of new mu.0
	mu.0 <- rbind(mu.0, S.tilde[idx,])
# sum(duplicated(mu.0))

  	###
	### Create receptacles for output
	###
  
	h.match.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
	h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
	theta.save <- numeric(n.mcmc)  # concentration parameter
	mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
	n.cls.save <- numeric(n.mcmc)
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	alpha.save <- matrix(0,n.mcmc,qW)
	beta.save <- matrix(0,n.mcmc,qX)
	gamma.save <- matrix(0,n.mcmc,qU)
	v.save <- matrix(0,T+n,n.mcmc)
	sigma.save <- numeric(n.mcmc)  # telemetry measurement error
	sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	keep <- list(sigma=0,sigma.mu=0,mu.0=0,gamma=0)
    
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) cat(k,"");flush.console()

		###
		### Updates pertaining to wet/dry status of y and z
		###

		#  Sample v(t.tilde) (auxilliary variable for y) 
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# Sample v(t) (auxilliary variable for s) 
		z1 <- which(z==1)
		z0 <- which(z==0)
		z1.sum <- length(z1)
		z0.sum <- length(z0)
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,z1.sum)
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,z0.sum)
		
		#  Sample alpha (coefficients on basis expansion W)
		A.inv <- solve(t(W.comb)%*%W.comb+Sigma.alpha.inv)
		b <- t(W.comb)%*%(v-X.comb%*%beta)
		alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
	
		#  Sample beta (coefficients on covariates influencing P(y==1))
		A.inv <- solve(t(X.comb)%*%X.comb+Sigma.beta.inv)
	  	b <- t(X.comb)%*%(v-W.comb%*%alpha)  # +mu.beta%*%Sigma.beta.inv
	  	beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	    # Sample z(t) (prediction of haul-out indicator variable for times t)
		linpred <- X.comb%*%beta+W.comb%*%alpha
		p <- pnorm(linpred[1:T,])
		p.tmp1 <- p*dnorm(s[,1],mu.0[h.match,1],sigma,log=FALSE)*
			dnorm(s[,2],mu.0[h.match,2],sigma,log=FALSE)
		p.tmp2 <- (1-p)*dnorm(s[,1],mu.0[h.match,1],sd.tmp,log=FALSE)*
			dnorm(s[,2],mu.0[h.match,2],sd.tmp,log=FALSE)
		p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
		z <- rbinom(T,1,p.tmp)
# z <- start$z
		
		###
	    ### Sample mu.0 (true location of clusters)
		### Note: sampling order does not matter here
	    ###
# browser()			
		# Sample 'occupied' mu.0 (clusters with non-zero membership)
				   
		# Use for data.table functionality
		# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0
			# get.mu.0(x,dt.h.idx[1:T,h.idx],z,s,mu,sigma,sigma.mu,S.tilde)))

		# Use for base functionality	
		mu.0.star <- matrix(NA,n.cls,2)  # ordered same as cls.idx
		S.match.star <- S.match+sample(tune$mu.0,n.cls,replace=TRUE) # ordered same as cls.idx
		idx <- which(S.match.star>0&S.match.star<max(S.tilde.idx)&
			!duplicated(c(S.match,S.match.star))[-(1:n.cls)])  # unique mu.0.star in S.tilde
		if(length(idx)>0){
			mu.0.star[idx,] <- S.tilde[S.match.star[idx],]
			mh.mu.0 <- t(sapply(idx,function(x)  # proposals for mu.0	
				get.mh.mu.0(x,s,z,h.match,cls.idx,mu.0[cls.idx,],mu.0.star,sigma,sigma.mu,
				S.match,S.match.star,U,gamma)))   
			idx <- idx[mh.mu.0]
			keep$mu.0 <- keep$mu.0+length(idx)
			mu.0[cls.idx[idx],] <- mu.0.star[idx,]  # accept proposals in S.tilde
			S.match[idx] <- S.match.star[idx]
		}
		
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|gamma]
# if(k==1000) browser()
		idx <- S.tilde.idx[-S.match]
		p <- exp(U[idx,]%*%gamma)
		# p <- exp(U[idx,-1]*gamma[-1])
		idx <- sample(idx,H-n.cls,replace=FALSE,prob=p)  # idx of new mu.0
# track proposals for mu.0
# points(S.tilde[idx,],pch=19,col=4,cex=0.25)
		mu.0[-cls.idx,] <- S.tilde[idx,]


		###
	    ### Sample gamma (selection coefficients for haul-out sites mu_0)
	    ###
# browser()
		gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		gamma.int.star <- log(sum(exp(U%*%gamma.star)))
		gamma.int <- log(sum(exp(U%*%gamma)))
		mh.star.gamma <- sum(U[S.match,]%*%gamma.star) - n.cls*gamma.int.star +
			sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))
		mh.0.gamma <- sum(U[S.match,]%*%gamma) - n.cls*gamma.int +
			sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))
		if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    gamma <- gamma.star
	        gamma.int <- gamma.int.star
	        keep$gamma <- keep$gamma+1
    	} 


		###
	    ### Sample sigma (observation error)
	    ###
# browser()
	    sigma.star <- rnorm(1,sigma,tune$sigma)
	    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
			idx <- z==1
			mu.0.tmp <- mu.0[h.match,]
	    	sd.tmp.star <- sqrt(sigma.star^2+sigma.mu^2)
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
				sd.tmp <- sd.tmp.star
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
			mu.0.tmp <- mu.0[h.match[idx],]
	    	sd.tmp.star <- sqrt(sigma^2+sigma.mu.star^2)
		    mh.star.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sd.tmp.star,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sd.tmp.star,log=TRUE))
		    mh.0.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sd.tmp,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sd.tmp,log=TRUE))
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
				sd.tmp <- sd.tmp.star
	        	Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		        keep$sigma.mu <- keep$sigma.mu+1
	    	} 
	    }
# sigma.mu <- start$sigma.mu

	
		###
		### Dirichlet process parameters
		### Note: sampling order matters here. Cluster parameters must be 
		### sampled in same order as pie, i.e., sorted by decreasing membership
		###
# browser()	
		# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# and Gelman et al. (2014), Section 23.3
		
	    # Sample h.t (cluster assignment indicator)
		mu.0.tmp <- mu.0[cls.samp,]  # same order as pie, i.e., by decreasing membership
		sd.tmp <- sqrt(sigma^2+sigma.mu^2)
	   	h.match <- sapply(1:T,function(x) sample(cls.samp,1,
			prob=pie*(dnorm(s[x,1],mu.0.tmp[,1],sigma)*
			dnorm(s[x,2],mu.0.tmp[,2],sigma))^z[x]*
			(dnorm(s[x,1],mu.0.tmp[,1],sd.tmp)*
			dnorm(s[x,2],mu.0.tmp[,2],sd.tmp))^(1-z[x])))

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
		cls.tab <- table(h.match)  # tabulate cluster membership
		n.cls <- length(cls.tab)  # number of clusters
		cls.ord <- order(cls.tab,decreasing=TRUE)  # clusters ordered by membership
		cls.idx <- as.numeric(names(cls.tab))  # idx of occupied clusters
		cls.diff <- setdiff(1:H,cls.idx)  # idx of unoccupied clusters
		cls.samp <- c(cls.idx[cls.ord],cls.diff)  # order in which clusters are sampled
		S.match <- S.tilde.idx[match(mu.0[cls.idx,2],S.tilde[,2])]  # S.tilde membership 
			# indicator

 	    # Stick-breaking process
    
	    # Use for data.table functionality
		# tab.cls.tmp <- c(rev(dt.tab.cls[,N]),rep(0,H-n.cls-1)) 
		
		# Use for base functionality
 	  	cls.tab.tmp <- c(cls.tab[cls.ord],rep(0,H-n.cls-1))  # membership in 
 	  		# decreasing order
 	  	# tab.cls.tmp <- c(tab.cls,rep(0,H-n.cls-1))  # membership in decreasing order
		eta <- c(rbeta(H-1,1+cls.tab.tmp,theta+T-cumsum(cls.tab.tmp)),1)  # stick-breaking
	    	# weights
	    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # Sample theta (concentration parameter); See Gelman section 23.3
       	theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta

			
		###
		###  Save samples 
		###
		
		# Use with data.table functionality
		# h.idx.save[,k] <- dt.h.idx[1:T,h.idx]
		# h.save[,,k] <- mu.0[dt.h.idx[1:T,h.idx],]
		
		# Use with base functionality
		h.match.save[,k] <- h.match 
		h.save[,,k] <- mu.0[h.match,]
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
	keep$gamma <- keep$gamma/n.mcmc
	cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	cat(paste("\ngamma acceptance rate:",round(keep$gamma,2))) 
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	list(h.match=h.match.save,h=h.save,mu.0=mu.0.save,beta=beta.save,alpha=alpha.save,
		gamma=gamma.save,theta=theta.save,sigma=sigma.save,sigma.mu=sigma.mu.save,z=z.save,
		v=v.save,n.cls=n.cls.save,keep=keep,n.mcmc=n.mcmc)
}