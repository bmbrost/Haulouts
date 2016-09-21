haulouts.1.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S.tilde,sigma.alpha,
	priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (04 SEP 2015)
 	### See haulouts.1.sim.R to simulate data according to this model specification,
 	### and haulouts.1.pdf for the model description, model statement, and
 	### full conditional distributions
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s; X.tilde=design matrix containing covariates
 	### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	### W.tilde=basis expansion for y; S.tilde=support Dirichlet process mixture 
 	### (i.e., haul-out sites); sigma.alpha=standard deviation of parameter 
 	### model for 'random' effects...
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
  
	get.mh.mu.0 <- function(x,mu.0,mu.0.star,S.tilde,ht,s,z,sigma,sigma.z0){
		# browser()
		mu.0.tmp <- S.tilde[mu.0[x],3:4]
		mu.0.tmp.star <- S.tilde[mu.0.star[x],3:4]
		idx.0 <- which(ht==mu.0[x]&z==0)
		idx.1 <- which(ht==mu.0[x]&z==1)
		# U.0 <- U[ht[x],]
		# U.star <- U[ht.star[x],]
		mh.star <- sum(dnorm(s[idx.1,1],mu.0.tmp.star[1],sigma,log=TRUE)+
			dnorm(s[idx.1,2],mu.0.tmp.star[2],sigma,log=TRUE),
			dnorm(s[idx.0,1],mu.0.tmp.star[1],sigma.z0,log=TRUE)+
			dnorm(s[idx.0,2],mu.0.tmp.star[2],sigma.z0,log=TRUE))
		mh.0 <- sum(dnorm(s[idx.1,1],mu.0.tmp[1],sigma,log=TRUE)+
			dnorm(s[idx.1,2],mu.0.tmp[2],sigma,log=TRUE),
			dnorm(s[idx.0,1],mu.0.tmp[1],sigma.z0,log=TRUE)+
			dnorm(s[idx.0,2],mu.0.tmp[2],sigma.z0,log=TRUE))
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
	###  Setup Variables 
	###
  
# browser() 
	T <- nrow(s)  # number of telemetry locations
	n <- length(y)  # number of wet/dry observations
	qX <- ncol(X)
	qW <- ncol(W)
	v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	X.comb <- rbind(X,X.tilde)  # combined design matrix for updates on alpha, beta
	W.comb <- rbind(W,W.tilde)  # combined design matrix for updates on alpha, beta 
	idx <- which(values(S.tilde)==1)
	S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))


	###
	### Starting values and priors
	###
# browser() 
	beta <- matrix(start$beta,qX)
	alpha <- matrix(0,qW)
	# alpha <- start$alpha
	theta <- start$theta
	sigma <- start$sigma
	sigma.mu <- start$sigma.mu
	pie <- start$pie
	z <- start$z
	H <- priors$H
	
	Sigma.inv <- solve(sigma^2*diag(2))
	Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
  	mu.beta <- matrix(0,qX,1)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	mu.alpha <- matrix(0,qW,1)
	Sigma.alpha <- diag(qW)*sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)
	sigma.z0 <- sqrt(sigma^2+sigma.mu^2)

	y1 <- which(y==1)+T
	y0 <- which(y==0)+T
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	linpred <- X.comb%*%beta+W.comb%*%alpha

	
	###
	### Set up Dirichlet process mixture variables
	###

	ht <- match(start$ht,S.tilde[,2])

	# Tabulate cluster membership with base functions
	tab <- table(ht)  # tabulate cluster membership
	m <- length(tab)  # number of clusters
	ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	mu.0 <- as.numeric(names(tab))  # idx of occupied clusters

  	
  	###
	### Create receptacles for output
	###
  
	ht.save <- matrix(0,T,n.mcmc)
	theta.save <- numeric(n.mcmc)  # concentration parameter
	m.save <- numeric(n.mcmc)
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	v.save <- matrix(0,T+n,n.mcmc)
	alpha.save <- matrix(0,n.mcmc,qW)
	beta.save <- matrix(0,n.mcmc,qX)
	sigma.save <- numeric(n.mcmc)  # telemetry measurement error
	sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	keep <- list(mu.0=0,sigma=0,sigma.mu=0,theta=0)
    
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) {
	    	cat(k,"");flush.console()	
    		plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	} 

		###
		### Dirichlet process parameters
		### Note: sampling order matters here. Cluster parameters must be 
		### sampled in same order as pie, i.e., sorted by decreasing membership
		###
	
		# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# and Gelman et al. (2014), Section 23.3

		# Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		   
		# Note: sampling order does not matter here

# browser()			
		# plot(S.tilde[,3:4],pch=19,cex=0.5,col="grey85")
		# points(S.tilde[mu.0,3:4],pch=19,col=1)
		mu.0.star <- sapply(mu.0,function(x)  # proposals for mu.0	
			sample(S.tilde[,1],1,prob=dnorm(S.tilde[x,3],S.tilde[,3],tune$mu.0)*
			dnorm(S.tilde[x,4],S.tilde[,4],tune$mu.0)))  
		# points(S.tilde[mu.0.star,3:4],col=rgb(1,0,0,0.25),cex=0.5,pch=19)
		dup.idx <- which(!duplicated(mu.0.star))  # exclude duplicate proposals

		mh.mu.0 <- sapply(dup.idx,function(x)  # accepted proposals
			get.mh.mu.0(x,mu.0,mu.0.star,S.tilde,ht,s,z,sigma,sigma.mu)) 
		keep$mu.0 <- keep$mu.0+sum(mh.mu.0)
		keep.idx <- dup.idx[mh.mu.0]
		mu.0[keep.idx] <- mu.0.star[keep.idx]
			
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|S.tilde]
		idx <- sample(S.tilde[-mu.0,1],H-m,replace=FALSE)  # idx of new mu.0
		# sum(idx%in%ht)
	    
	    # Sample cluster assignment indicator, ht (Note: sampling order matters here)
		samp <- c(mu.0[ord],idx)  # sampling in order of decreasing membership
		# sum(duplicated(samp))
	   	# ht <- sapply(1:T,function(x) sample(samp,1,
			# prob=pie*(dnorm(s[x,1],S.tilde[samp,3],sigma)*
			# dnorm(s[x,2],S.tilde[samp,4],sigma))^z[x]*
			# (dnorm(s[x,1],S.tilde[samp,3],sigma.z0)*
			# dnorm(s[x,2],S.tilde[samp,4],sigma.z0))^(1-z[x])))
		ht <- sapply(1:T,function(x) sample(samp,1,prob= 
			exp(log(pie)+z[x]*(dnorm(s[x,1],S.tilde[samp,3],sigma,log=TRUE)+
			dnorm(s[x,2],S.tilde[samp,4],sigma,log=TRUE))+	
			(1-z[x])*(dnorm(s[x,1],S.tilde[samp,3],sigma.z0,log=TRUE)+
			dnorm(s[x,2],S.tilde[samp,4],sigma.z0,log=TRUE)))))
		
		# Tabulate cluster membership with base functions
		tab <- table(ht)  # tabulate cluster membership
		m <- length(tab)  # number of clusters
		ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
		mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		
	    # Stick-breaking process (Note: sampling order matters here)
		tab.tmp <- c(tab[ord],rep(0,H-m-1))  # membership in decreasing order
		eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # Sample theta (concentration parameter); See Gelman section 23.3, 
	    # Ishwaran and Zarepour (2000)
       	# theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta

		# Sample theta (concentration parameter); See Escobar and West (1995) and 
		# West (1997?) white paper on hyperparameter estimation in DP
		tmp <- rbeta(1,theta+1,T)
		c <- priors$r
		d <- priors$q
		p.tmp <- (c+m-1)/(c+m-1+T*(d-log(tmp)))
		p.tmp <- rbinom(1,1,p.tmp)
		theta <- ifelse(p.tmp==1,rgamma(1,c+m,d-log(tmp)),rgamma(1,c+m-1,d-log(tmp)))

		# Metropolis-Hastings update based on uniform prior on theta
		# theta.star <- rnorm(1,theta,0.5)
	    # if(theta.star>0 & theta.star<10){
	    	# mh.star.theta <- m*log(theta.star)+lgamma(theta.star)-lgamma(theta.star+T)
	    	# mh.0.theta <- m*log(theta)+lgamma(theta)-lgamma(theta+T)
		    # if(exp(mh.star.theta-mh.0.theta)>runif(1)){
	    	    # theta <- theta.star
		        # keep$theta <- keep$theta+1
	    	# } 
	    # }

	
		###
		### Updates pertaining to wet/dry status of y and z
		###
# browser()
		# Sample v(t.tilde) (auxilliary variable for y) 
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# Sample v(t) (auxilliary variable for s) 
		z1 <- which(z==1)
		z0 <- which(z==0)
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))

		#  Sample alpha (coefficients on basis expansion W)
		A.inv <- solve(t(W.comb)%*%W.comb+Sigma.alpha.inv)
		b <- t(W.comb)%*%(v-X.comb%*%beta)
		alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
		
		#  Sample beta (coefficients on covariates influencing haul-out probability)
		A.inv <- solve(t(X.comb)%*%X.comb+Sigma.beta.inv)
	  	b <- t(X.comb)%*%(v-W.comb%*%alpha)  # +mu.beta%*%Sigma.beta.inv
	  	beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	    # Sample z(t) (prediction of haul-out indicator variable for times t)
		linpred <- X.comb%*%beta+W.comb%*%alpha
		p <- pnorm(linpred[1:T,])
		p.tmp1 <- p*dnorm(s[,1],S.tilde[ht,3],sigma,log=FALSE)*
			dnorm(s[,2],S.tilde[ht,4],sigma,log=FALSE)
		p.tmp2 <- (1-p)*dnorm(s[,1],S.tilde[ht,3],sigma.z0,log=FALSE)*
			dnorm(s[,2],S.tilde[ht,4],sigma.z0,log=FALSE)
		# p.tmp <- p.tmp1/(p.tmp1+p.tmp2)
		p.tmp <- exp(log(p.tmp1)-log(p.tmp1+p.tmp2))	
		z <- rbinom(T,1,p.tmp)
# z <- start$z

	
	    ###
	    ### Sample sigma (observation error)
	    ###
# browser()
	    sigma.star <- rnorm(1,sigma,tune$sigma)
	    if(sigma.star>priors$sigma.l & sigma.star<priors$sigma.u){
			idx <- z==1
			mu.0.tmp <- S.tilde[ht,3:4]
	    	sigma.z0.star <- sqrt(sigma.star^2+sigma.mu^2)
	    	mh.star.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma.star,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma.star,log=TRUE))+
		    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sigma.z0.star,log=TRUE)+
		    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sigma.z0.star,log=TRUE))
		    mh.0.sigma <- sum(dnorm(s[idx,1],mu.0.tmp[idx,1],sigma,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[idx,2],sigma,log=TRUE))+
		    	sum(dnorm(s[!idx,1],mu.0.tmp[!idx,1],sigma.z0,log=TRUE)+
		    	dnorm(s[!idx,2],mu.0.tmp[!idx,2],sigma.z0,log=TRUE))
		    if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
	    	    sigma <- sigma.star
				sigma.z0 <- sigma.z0.star
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
	    if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			idx <- which(z==0)
			mu.0.tmp <- S.tilde[ht[idx],3:4]
	    	sigma.z0.star <- sqrt(sigma^2+sigma.mu.star^2)
		    mh.star.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sigma.z0.star,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sigma.z0.star,log=TRUE))
		    mh.0.sigma.mu <- sum(dnorm(s[idx,1],mu.0.tmp[,1],sigma.z0,log=TRUE)+
		    	dnorm(s[idx,2],mu.0.tmp[,2],sigma.z0,log=TRUE))
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
        		sigma.z0 <- sigma.z0.star
	        	Sigma.mu.inv <- solve(sigma.mu^2*diag(2))
		        keep$sigma.mu <- keep$sigma.mu+1
	    	} 
	    }
# sigma.mu <- start$sigma.mu
	
		###
		###  Save samples 
		###
		
# browser()	
		ht.save[,k] <- S.tilde[ht,2]
		theta.save[k] <- theta    
		sigma.save[k] <- sigma
		sigma.mu.save[k] <- sigma.mu
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		v.save[,k] <- v
		z.save[,k] <- z
		m.save[k] <- m
	}
  
	###
	### Write output
	###
	  
	keep$sigma <- keep$sigma/n.mcmc
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu.0 <- keep$mu.0/sum(m.save)
	cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
cat(paste("\ntheta acceptance rate:",round(keep$theta/n.mcmc,2))) 

	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	list(ht=ht.save,beta=beta.save,alpha=alpha.save,
		theta=theta.save,sigma=sigma.save,sigma.mu=sigma.mu.save,z=z.save,v=v.save,
	  	m=m.save,keep=keep,n.mcmc=n.mcmc)
}