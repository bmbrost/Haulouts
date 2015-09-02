haulouts.1.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S.tilde,sigma.alpha,
	priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (13 AUG 2015)
 	### See haulouts.sim.R to simulate data according to this model specification,
 	### and haulouts.pdf for the model description, model statement, and
 	### full conditional distributions
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s; X.tilde=design matrix containing covariates
 	### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	### W.tilde=basis expansion for y; S.tilde=support Dirichlet process mixture 
 	### (i.e., haul-out sites); sigma.alpha=standard deviation of parameter 
 	### model for 'random' effects
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
	v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	X.comb <- rbind(X,X.tilde)  # combined design matrix for updates on alpha, beta
	W.comb <- rbind(W,W.tilde)  # combined design matrix for updates on alpha, beta 
	idx <- which(values(S.tilde)>0)
	S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))

	# S.tilde.idx <- which(!is.na(values(S.tilde)))
	# S.tilde.xy <- xyFromCell(S.tilde,S.tilde.idx)
	# S.tilde.idx <- values(S.tilde)[S.tilde.idx]


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
	sd.tmp <- sqrt(sigma^2+sigma.mu^2)

	y1 <- which(y==1)+T
	y0 <- which(y==0)+T
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	linpred <- X.comb%*%beta+W.comb%*%alpha

	###
	### Set up Dirichlet process mixture variables
	###

	# mu.0 <- unique(h)  # unique cluster location s
	# n.cls <- nrow(mu.0)  # number of clusters
	# h.match <- match(paste(start$h[,1],start$h[,2]),
		# paste(mu.0[,1],mu.0[,2]))  # cluster membership indicator

	S.match <- match(paste(start$h[,1],start$h[,2]),
		paste(S.tilde[,3],S.tilde[,4]))  # S.tilde membership indicator


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
	# cls.tab <- table(h.match)  # tabulate cluster membership
	# cls.tab <- table(S.match)  # tabulate cluster membership
	# n.cls <- length(cls.tab)  # number of clusters
	# # cls.ord <- order(cls.tab,decreasing=TRUE)  # clusters ordered by membership
	# cls.idx <- as.numeric(names(cls.tab))  # idx of occupied clusters
	
	cls.idx <- unique(S.match)
	n.cls <- length(cls.idx)
	
	# cls.diff <- setdiff(1:H,cls.idx)  # idx of unoccupied clusters
	# cls.samp <- c(cls.idx[cls.ord],cls.diff)  # order in which clusters are sampled
		
	# S.match <- match(paste(mu.0[cls.idx,1],mu.0[cls.idx,2]),
		# paste(S.tilde.xy[,1],S.tilde.xy[,2]))  # S.tilde membership indicator
  
	# Propose values for mu.0
	# idx <- sample(S.tilde.idx[-S.match],H-n.cls,replace=FALSE)  # idx of new mu.0
	# mu.0 <- rbind(mu.0, S.tilde.xy[idx,])
# sum(duplicated(mu.0))


  	###
	### Create receptacles for output
	###
  
	# h.match.save <- matrix(0,T,n.mcmc)  # cluster assignment indicator variable
	# h.save <- array(0,dim=c(T,2,n.mcmc))  # cluster assignment indicator variable
	S.match.save <- matrix(0,T,n.mcmc)
	theta.save <- numeric(n.mcmc)  # concentration parameter
	# mu.0.save <- array(0,dim=c(H,2,n.mcmc))  # cluster locations
	n.cls.save <- numeric(n.mcmc)
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	alpha.save <- matrix(0,n.mcmc,qW)
	beta.save <- matrix(0,n.mcmc,qX)
	v.save <- matrix(0,T+n,n.mcmc)
	sigma.save <- numeric(n.mcmc)  # telemetry measurement error
	sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	keep <- list(mu.0=0,sigma=0,sigma.mu=0)
    
	###
	### Begin MCMC loop
	###
  
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) {
	    	cat(k,"");flush.console()	
    		plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	} 

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
			z1.sum <- length(z1)
			z0.sum <- length(z0)
		  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,z1.sum)
			v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,z0.sum)
	
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
			p.tmp1 <- p*dnorm(s[,1],mu.0[h.match,1],sigma,log=FALSE)*
				dnorm(s[,2],mu.0[h.match,2],sigma,log=FALSE)
			p.tmp2 <- (1-p)*dnorm(s[,1],mu.0[h.match,1],sd.tmp,log=FALSE)*
				dnorm(s[,2],mu.0[h.match,2],sd.tmp,log=FALSE)
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

		###
	    ### Sample mu.0 (true location of occupied clusters)
		### Note: sampling order does not matter here
	    ###

		# Sample 'occupied' mu.0 (clusters with non-zero membership)		   
	
		# Use for data.table functionality
		# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0
			# get.mu.0(x,dt.h.idx[1:T,h.idx],z,s,mu,sigma,sigma.mu,S.tilde)))
browser()			

plot(S.tilde)
points(mu.0[cls.idx,],pch=19)
		# Proposals for mu.0 sampled from S.tilde
		# S.match.star <- sapply(cls.idx,function(x)  # ordered same as cls.idx
			# sample(S.tilde.idx,1,prob=dnorm(mu.0[x,1],S.tilde.xy[,1],tune$mu.0)*
			# dnorm(mu.0[x,2],S.tilde.xy[,2],tune$mu.0)))
		# mu.0.star <- S.tilde.xy[S.match.star,]  # ordered same as cls.idx
		# points(mu.0.star,col=rgb(1,0,0,0.25),cex=0.5,pch=19)

		cls.idx.star <- sapply(cls.idx,function(x)  # ordered same as cls.idx
			sample(S.tilde[,1],1,prob=dnorm(S.tilde[x,3],S.tilde[,3],tune$mu.0)*
			dnorm(S.tilde[x,4],S.tilde[,4],tune$mu.0)))
		# points(S.tilde.xy[S.match.star,],col=rgb(1,0,0,0.25),cex=0.5,pch=19)

		# mu.0.star <- S.tilde.xy[S.match.star,]  # ordered same as cls.idx
# points(mu.0.star,col=rgb(1,0,0,0.25),cex=0.5,pch=19)

	get.mh.mu.0 <- function(x,cls.idx,cls.idx.star,S.tilde,S.match,s,z,sigma,sigma.mu){
		# ,U,gamma){
		# browser()
		# x
		mu.0 <- S.tilde[cls.idx[x],3:4]
		mu.0.star <- S.tilde[cls.idx.star[x],3:4]
		idx.0 <- which(S.match==cls.idx[x]&z==0)
		idx.1 <- which(S.match==cls.idx[x]&z==1)
    	sd.tmp <- sqrt(sigma^2+sigma.mu^2)
		# U.0 <- U[S.match[x],]
		# U.star <- U[S.match.star[x],]
		mh.star <- sum(dnorm(s[idx.1,1],mu.0.star[1],sigma,log=TRUE)+
			dnorm(s[idx.1,2],mu.0.star[2],sigma,log=TRUE))+
			sum(dnorm(s[idx.0,1],mu.0.star[1],sd.tmp,log=TRUE)+
			dnorm(s[idx.0,2],mu.0.star[2],sd.tmp,log=TRUE))  # +			
			# U.star%*%gamma # - int
		mh.0 <- sum(dnorm(s[idx.1,1],mu.0[1],sigma,log=TRUE)+
			dnorm(s[idx.1,2],mu.0[2],sigma,log=TRUE))+
			sum(dnorm(s[idx.0,1],mu.0[1],sd.tmp,log=TRUE)+
			dnorm(s[idx.0,2],mu.0[2],sd.tmp,log=TRUE))  # +			
			# U.0%*%gamma # - int
		exp(mh.star-mh.0)>runif(1)
	}

cls.idx
S.match.star

		mh.mu.0 <- sapply(1:n.cls,function(x)  # proposals for mu.0	
			get.mh.mu.0(x,cls.idx,cls.idx.star,S.tilde,S.match,s,z,sigma,sigma.mu)) 
			#,U,gamma)))   
		keep$mu.0 <- keep$mu.0+sum(mh.mu.0)
		mu.0[cls.idx[mh.mu.0],] <- mu.0.star[mh.mu.0,]  # accept proposals in S.tilde
		S.match[mh.mu.0] <- S.match.star[mh.mu.0]






	# get.mh.mu.0 <- function(x,s,z,S.match,cls.idx,mu.0,mu.0.star,sigma,sigma.mu,
		# S.match,S.match.star){  #,U,gamma){
		# # browser()
		# # x
		# # mu.0.tmp <- mu.0[x,]
		# # mu.0.star.tmp <- mu.0.star[x,]
		# idx.0 <- which(h.match==cls.idx[x]&z==0)
		# idx.1 <- which(h.match==cls.idx[x]&z==1)
    	# sd.tmp <- sqrt(sigma^2+sigma.mu^2)
		# # U.0 <- U[S.match[x],]
		# # U.star <- U[S.match.star[x],]
		# mh.star <- sum(dnorm(s[idx.1,1],mu.0.star[x,1],sigma,log=TRUE)+
			# dnorm(s[idx.1,2],mu.0.star[x,2],sigma,log=TRUE))+
			# sum(dnorm(s[idx.0,1],mu.0.star[x,1],sd.tmp,log=TRUE)+
			# dnorm(s[idx.0,2],mu.0.star[x,2],sd.tmp,log=TRUE))  # +			
			# # U.star%*%gamma # - int
		# mh.0 <- sum(dnorm(s[idx.1,1],mu.0[x,1],sigma,log=TRUE)+
			# dnorm(s[idx.1,2],mu.0[x,2],sigma,log=TRUE))+
			# sum(dnorm(s[idx.0,1],mu.0[x,1],sd.tmp,log=TRUE)+
			# dnorm(s[idx.0,2],mu.0[x,2],sd.tmp,log=TRUE))  # +			
			# # U.0%*%gamma # - int
		# exp(mh.star-mh.0)>runif(1)
	# }

		# mh.mu.0 <- sapply(1:n.cls,function(x)  # proposals for mu.0	
			# get.mh.mu.0(x,s,z,h.match,cls.idx,mu.0[cls.idx,],mu.0.star,sigma,sigma.mu,
			# S.match,S.match.star))  #,U,gamma)))   
		# keep$mu.0 <- keep$mu.0+sum(mh.mu.0)
		# mu.0[cls.idx[mh.mu.0],] <- mu.0.star[mh.mu.0,]  # accept proposals in S.tilde
		# S.match[mh.mu.0] <- S.match.star[mh.mu.0]


# Use for base functionality	
# mu.0.tmp <- t(sapply(idx.cls,function(x)  # proposals for mu.0	
	# get.mu.0(x,h.match,z,s,Sigma.inv,Sigma.mu.inv)))  
# idx <- which(mu.0.tmp[,1]>S.tilde[1,1]&mu.0.tmp[,1]<S.tilde[2,1]&
	# mu.0.tmp[,2]>S.tilde[1,2]&mu.0.tmp[,2]<S.tilde[3,2])  # idx of mu.0 in S.tilde	
# mu.0[idx.cls[idx],] <- mu.0.tmp[idx,]  # accept proposals in S.tilde
# n.cls.star <- H-n.cls  # number of new clusters to propose
# mu.0[-idx.cls,] <- cbind(runif(n.cls.star,S.tilde[1,1],S.tilde[2,1]),
      # runif(n.cls.star,S.tilde[1,2],S.tilde[3,2]))  # update mu.0 with mu.star

		
		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|gamma]
# if(k==1000) browser()
		idx <- sample(S.tilde.idx[-S.match],H-n.cls,replace=FALSE)  # idx of new mu.0
		mu.0[-cls.idx,] <- S.tilde.xy[idx,]
# sum(duplicated(mu.0))
















			
# browser()
		    # Sample S.match (cluster assignment indicator)

			# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|gamma]
# if(k==1000) browser()
cls.idx
			cls.samp <- cls.idx[cls.ord]
		
			idx <- sample(S.tilde[-cls.idx,1],H-n.cls,replace=FALSE)  # idx of new mu.0


sum(idx%in%S.match)
			cls.samp <- c(cls.samp,idx)

			# mu.0.tmp <- rbind(S.tilde.xy[cls.idx[cls.ord],],S.tilde.xy[idx,])

# sum(duplicated(mu.0))



		   	S.match <- sapply(1:T,function(x) sample(cls.samp,1,
				prob=pie*(dnorm(s[x,1],S.tilde[cls.samp,3],sigma)*
				dnorm(s[x,2],S.tilde[cls.samp,4],sigma))^z[x]*
				(dnorm(s[x,1],S.tilde[cls.samp,3],sd.tmp)*
				dnorm(s[x,2],S.tilde[cls.samp,4],sd.tmp))^(1-z[x])))

		   	# h.match <- sapply(1:T,function(x) sample(cls.samp,1,
				# prob=pie*(dnorm(s[x,1],mu.0[cls.samp,1],sigma)*
				# dnorm(s[x,2],mu.0[cls.samp,2],sigma))^z[x]*
				# (dnorm(s[x,1],mu.0[cls.samp,1],sd.tmp)*
				# dnorm(s[x,2],mu.0[cls.samp,2],sd.tmp))^(1-z[x])))
		
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
			# cls.tab <- table(h.match)  # tabulate cluster membership
			cls.tab <- table(S.match)  # tabulate cluster membership
			n.cls <- length(cls.tab)  # number of clusters
			cls.ord <- order(cls.tab,decreasing=TRUE)  # clusters ordered by membership
			cls.idx <- as.numeric(names(cls.tab))  # idx of occupied clusters
			# cls.diff <- setdiff(1:H,cls.idx)  # idx of unoccupied clusters
			# cls.samp <- c(cls.idx[cls.ord],cls.diff)  # order in which clusters are sampled
			# S.match <- match(paste(mu.0[cls.idx,1],mu.0[cls.idx,2]),
				# paste(S.tilde.xy[,1],S.tilde.xy[,2]))  # S.tilde membership indicator
		  
 		    # Stick-breaking process
		    
		    # Use for data.table functionality
			# tab.cls.tmp <- c(rev(dt.tab.cls[,N]),rep(0,H-n.cls-1)) 
			
			# Use for base functionality
			cls.tab.tmp <- c(cls.tab[cls.ord],rep(0,H-n.cls-1))  # membership in 
 	  			# decreasing order
			eta <- c(rbeta(H-1,1+cls.tab.tmp,theta+T-cumsum(cls.tab.tmp)),1)  # stick-breaking
		    	# weights
		    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

		    # Sample theta (concentration parameter); See Gelman section 23.3
	       	theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta










	
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
	list(h.match=h.match.save,h=h.save,mu.0=mu.0.save,beta=beta.save,alpha=alpha.save,
		theta=theta.save,sigma=sigma.save,sigma.mu=sigma.mu.save,z=z.save,v=v.save,
	  	n.cls=n.cls.save,keep=keep,n.mcmc=n.mcmc)
}