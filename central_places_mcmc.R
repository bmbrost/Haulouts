# Software Disclaimer:
# Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, 
# expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning 
# of the program and related program material nor shall the fact of distribution constitute any 
# such warranty, and no responsibility is assumed by the USGS in connection therewith. 
#
#
# The MIT License (MIT)
#
# Copyright (c) 2016 Brian M. Brost
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#
# Function name: central.places.mcmc
#
# Description: Dirichlet process mixture model for estimating the location of central places from 
# satellite telemetry data.
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 03 OCT 2016
#
# Reference: Brost, B.M., M.B. Hooten, and R.J. Small. 2017. Leveraging constraints and 
# 	biotelemetry data to pinpoint repetitively used spatial features. Ecology.
#
# Required R packages: 
#
# Inputs:
# s  - matrix of telemetry locations s(t)
# lc - Argos location quality class corresponding to each row in s 
# y - ancillary behavioral data (i.e., binary wet/dry status, y(t)) corresponding to each row in s
# X - design matrix containing covariates influencing behavioral state y(t)
# W - basis expansion for modeling dependence in y(t)
# J - upper bound to the Dirichlet process truncation approximation
# S.tilde - matrix summarizing support of haul-out sites. Columns consist of:
#	1. index of rows in S.tilde (i.e., sequence of 1:nrow(S.tilde))
#	2. index of raster cells within S.tilde (i.e., cell numbers associated with S.tilde)
#	3. x-coordinate of raster cells within S.tilde
#	4. y-coordinate of raster cells within S.tilde
# priors - list of prior distribution parameters containing the following elements:
#	1. mu.sigma - mean of lognormal prior for sigma_mu
#	2. sigma.sigma - standard deviation of lognormal prior for sigma_mu
#	3. f.hat - base distribution of Dirichlet process
#	4. sigma.beta - standard deviation of normal prior for beta
#	5. u.sigma - upper limit of uniform prior for sigma
#	6. r.alpha - shape parameter of IG prior for sigma_alpha^2
#	7. q.alpha - rate parameter of IG prior for sigma_alpha^2
#	8. r.theta - shape parameter of Gamma prior for theta
#	9. q.theta - rate parameter of Gamma prior for theta
# tune - list of tuning parameters for Metropolis-Hastings updates containing the following elements:
#	1. mu - tuning parameter of proposal distribution for mu_j
#	2. sigma.mu - tuning parameter of proposal distribution for sigma_mu
#	3. sigma - vector of tuning parameters for sigma_c (one per Argos location quality class)
#	4. rho - vector of tuning parameters for rho_c (one per Argos location quality class)
# 	5. a - vector of tuning parameters for a_c (one per Argos location quality class)
# start - list of starting values containing the following elements:
#	1. mu - starting values for mu_j
#	2. pi - starting values for pi
#	3. theta - starting value for theta
#	4. sigma.mu - starting value for sigma_mu (standard deviation)
#	5. sigma - starting values for sigma_c (one per Argos location quality class)
#	6. rho - starting values for rho_c (one per Argos location quality class)
#	7. a - starting values for a_c (one per Argos location quality class)
#	8. alpha - starting values for alpha
#	9. sigma.alpha - starting value for sigma_alpha (standard deviation)
#	10. beta - starting values for beta
# adapt - enable adaptive tuning for Metropolis-Hastings updates? 
# n.mcmc - number of MCMC iterations desired
#
# Notes: The model statement and full-conditional distributions corresponding to 
#	this MCMC algorithm are presented in Appendix A of Brost et al. (2017). Associated 
#	pseudocode is presented in Appendix B.
#

central.places.mcmc <- function(s,lc,y,X,W,J,S.tilde,priors,tune,start,adapt=FALSE,n.mcmc){

	t.start <- Sys.time()
	cat(paste("\n\nStart time:",t.start,"\n"))


	#####################################################################
	### Libraries and Subroutines
	#####################################################################

	truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}
	
	tkern <- function(d,P,nu=100){  # kernel of t-distribution
		(1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	}

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}

	get.Sigma <- function(sigma2,a,rho){  # Get covariance matrix
		S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		list(P=P,b=b,S=S)
	}

	dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# Density of mixture t-distribution
		x <- matrix(x,,2)	
		y <- matrix(y,,2)
		if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		if(!is.null(lc)){
			S <- S[[lc]]
			Q <- Q[[lc]]
		}
		P <- ifelse(z==1,S,Q)[[1]]  # precision matrix
		b <- ifelse(z==1,S$b,Q$b)  # determinant
		d <- x-y
		out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
		if(!log) out <- 0.5*out*b^(-0.5)  # density
		out
	}
  
	get.mh.mu <- function(x,mu,mu.star,S,h,s,lc,z,Sigma,Q,f.hat){
		# Accept/reject proposals for mu
		mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		d <- log(f.hat[mu[x]])  # prior log density for mu
		d.star <- log(f.hat[mu.star[x]])  # prior log density for mu.star
		idx.0 <- which(h==x&z==0)  # obs. associated with mu and z=0
		idx.1 <- which(h==x&z==1)  # obs. associated with mu and z=1
		num.z1 <- denom.z1 <- num.z0 <- denom.z0 <- 0
		if(length(idx.1)>0){
			num.z1 <- sum(sapply(idx.1,function(x)
				dt2(s[x,],mu.xy.star,z=1,Sigma,Q,lc[x],log=TRUE)))
			denom.z1 <- sum(sapply(idx.1,function(x)
				dt2(s[x,],mu.xy,z=1,Sigma,Q,lc[x],log=TRUE)))
		}
		if(length(idx.0)>0){
			num.z0 <- sum(sapply(idx.0,function(x) 
				dt2(s[x,],mu.xy.star,z=0,Sigma,Q,lc[x],log=TRUE)))
			denom.z0 <- sum(sapply(idx.0,function(x)
				dt2(s[x,],mu.xy,z=0,Sigma,Q,lc[x],log=TRUE)))
		}
		mh.star <- num.z1+num.z0+d.star  # numerator of Metropolis ratio
		mh.0 <-	denom.z1+denom.z0+d  # denominator of Metropolis ratio
		exp(mh.star-mh.0)>runif(1)  # Accept or reject
	}


	#####################################################################
	###  Setup Variables 
	#####################################################################

	cat("\nSetting up variables....")
	
	Ts <- nrow(s)  # number of telemetry locations
	Ty <- length(y)  # number of wet/dry observations
	qX <- ncol(X)  # column dimension of X
	qW <- ncol(W)  # column dimension of W
	v <- numeric(Ty)  # auxilliary variable for haul-out process
	y1 <- which(y==1)  # idx of dry observations
	y0 <- which(y==0)  # idx of wet observations
	y1.sum <- length(y1)  # number of dry observations
	y0.sum <- length(y0)  # number of wet observations
	n.lc <- length(unique(lc))  # number of error classes
	lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
		# idx of observations by Argos location quality class
	y1.lc <- tapply(y1,lc[y1],I)  # by behavioral state and location quality class
	y0.lc <- tapply(y0,lc[y0],I)  # by behavioral state and location quality class
	lc <- as.numeric(lc)  # numerical Argos location quality class


	#####################################################################
	### Priors
	#####################################################################
	
	cat("\nGetting priors....")
	
	# Observation model
	u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma_mu
	sigma.sigma <- priors$sigma.sigma  # standard deviation of longormal prior for sigma_mu
	r.alpha <- priors$r.alpha  # IG prior on sigma_alpha^2
	q.alpha <- priors$q.alpha  # IG prior on sigma_alpha^2
	r.theta <- priors$r.theta  # Gamma prior on theta
	q.theta <- priors$q.theta  # Gamma prior on theta
	
	# Temporal process model 
	mu.alpha <- matrix(0,qW,1)  # mean of random effects
	mu.beta <- matrix(0,qX,1)  # mean of coefficients on X (fixed effects)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2  # variance for beta
	Sigma.beta.inv <- solve(Sigma.beta)
	A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# Spatial process model
	f.hat <- priors$f.hat  # prior on mu_j
	

	#####################################################################
	### Standardize parameters
	#####################################################################

	cat("\nStandardizing variables....")
	
	# Center and scale s and S.tilde
	s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		ncol=2,byrow=TRUE))/s.sd
	
	# Center and scale tuning parameters
	tune$sigma.mu <- tune$sigma.mu/s.sd
	tune$mu <- tune$mu/s.sd
	tune$sigma <- tune$sigma/s.sd

	# Center and scale priors
	mu.sigma <- mu.sigma/s.sd  
	u.sigma <- u.sigma/s.sd
	
	# Center and scale starting values
	sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	sigma <- start$sigma/s.sd  # observation model standard deviation


	#####################################################################
	### Starting values: Appendix B, Steps 1 and 2
	#####################################################################

	cat("\nGetting starting values....")

	# Observation model
	a <- start$a
	rho <- start$rho
	Sigma <- sapply(1:n.lc,function(x)
		get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	Q <- sapply(1:n.lc,function(x) 
		get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# Temporal process model
	alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	Sigma.alpha <- diag(qW)*start$sigma.alpha^2  # variance of random effects
  	Sigma.alpha.inv <- solve(Sigma.alpha)
	W.cross <- t(W)%*%W  # cross product of W
  	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	
	# Spatial process model
	theta <- start$theta  # DP concentration parameter
	mu.j <- match(start$mu,S.tilde[,2])
	m <- length(mu.j)  # number of occupied clusters
	mu.j <- c(mu.j,sample(S.tilde[-mu.j,1],J-m,replace=FALSE))
	pi <- c(start$pi,rep(0,J-m))  # probability mass for each mu
	mu.star <- numeric(J)  # receptacle for proposals


  	#####################################################################
	### Create receptacles for output
	#####################################################################
  
  	cat("\nCreating receptacles for output....")
	
	sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudinal telemetry measurement error
	a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
	beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients
	alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	theta.save <- numeric(n.mcmc)  # DP concentration parameter
	m.save <- numeric(n.mcmc)  # number of clusters
	v.save <- matrix(0,Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	sigma.mu.save <- numeric(n.mcmc)  # homerange dispersion parameter
	sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects
    
    
    #####################################################################
	### MCMC loop: Appendix B, Steps 3-9
	#####################################################################
	
	# Track overall MH accpetance rate
	keep <- list(mu=0,sigma.mu=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# Track MH accpetance rate for adaptive tuning
	keep.tmp <- keep  
	T.b <- 50  # frequency of adaptive tuning
	m.save.tmp <- 0  # number of clusters

	t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	t.mcmc.start <- Sys.time()  # timing MCMC iterations

	cat("\nEntering MCMC Loop....\n")
	
	for (k in 1:n.mcmc) {  # Appendix B, Step 9
    	
    	if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	cat(k,"");flush.console()	
    		plot(pi,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	} 
		
		if(adapt==TRUE & k%%50==0) {  # Adaptive tuning
			keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			tune$sigma.mu <- get.tune(tune$sigma.mu,keep.tmp$sigma.mu,k)
			tune$mu <- get.tune(tune$mu,keep.tmp$mu,k)
			tune$sigma <- sapply(1:n.lc,function(x) 
				get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			tune$a <- sapply(1:n.lc,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			tune$rho <- sapply(1:n.lc,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
			m.save.tmp <- 0
	   	} 	
	
		
		#--------------------------------------------------------------------------
		# Update spatial process model parameters: Appendix B, Step 4
		#--------------------------------------------------------------------------

		###
	    ### Update latent class status h(t): Appendix B, Step 4(a)
	    ###

		# Indexes records in mu.j (1:J)
		h <- sapply(1:Ts,function(x) sample(1:J,1,prob= 
			exp(log(pi)+dt2(s[x,],S.tilde[mu.j,3:4],y[x],Sigma,Q,lc[x]))))

		###
		### Tabulate cluster membership m_j: Appendix B, Step 4(b) 
		###

		n <- table(h)
		m <- length(n)  # number of clusters			
		j <- as.numeric(names(n))  # idx of occupied clusters
	
		###
		### Update the stick-breaking weights eta_j: Appendix B, Step 4(c)
		###

		n.tmp <- numeric(J)
		n.tmp[j] <- n
		eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights
		eta[-J] <- abs(eta[-J]-10^-5)  # fudge factor for preventing eta[1:(J-1)]=1 or 0

		###
		### Update the stick-breaking probabilities pi_j: Appendix B, Step 4(d)
		###
		
	    pi <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities
	
		###
		### Update the DP concentration parameter theta: Appendix B, Step 4(e)
		###

		theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
		
		###
		### Update 'occupied' mu_j: Appendix B, Step 4(f)		   
		###

		# mu.j indexes rows in S.tilde
		mu.star[j] <- sapply(mu.j[j],
			function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*
			dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		idx <- j[which(!duplicated(mu.star[j]))]  # exclude duplicate proposals
		mh <- sapply(idx,function(x)  # accepted proposals
			get.mh.mu(x,mu.j,mu.star,S.tilde,h,s,lc,y,Sigma,Q,f.hat)) 
		keep$mu <- keep$mu+sum(mh)
		keep.tmp$mu <- keep.tmp$mu+sum(mh)
		keep.idx <- idx[mh]
		mu.j[keep.idx] <- mu.star[keep.idx]  # update mu

		###
		### Update 'unoccupied' mu_j: Appendix A, Step 4(g)
		###

		mu.j[-j] <- sample(S.tilde[-mu.j[j],1],J-m,replace=FALSE,prob=f.hat[-mu.j[j]])  
			# idx of unoccupied mu

		###
		### Use h to map mu_j to s(t): Appendix B, Step 4(h)
		###

		mu.t <- S.tilde[mu.j[h],3:4]  # row idx of S.tilde


		#--------------------------------------------------------------------------
	  	# Update movement parameter sigma_mu: Appendix B, Step 5 
	    #--------------------------------------------------------------------------

	    sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		if(sigma.mu.star>0){
			Q.star <- sapply(1:n.lc,function(x) 
				get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
		    mh.star.sigma.mu <-	sum(sapply(y0,function(x) 
		    	dt2(s[x,],mu.t[x,],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    mh.0.sigma.mu <- sum(sapply(y0,function(x)
		    	dt2(s[x,],mu.t[x,],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
				Q <- Q.star
				keep$sigma.mu <- keep$sigma.mu+1
				keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	} 
		}
	
	
		#--------------------------------------------------------------------------
		# Update observation model parameters: Appendix B, Step 6
		#--------------------------------------------------------------------------

		sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		for(i in 1:n.lc){  # iterate over error classes: Appendix B, Step 6(e)

			# Subset data for Argos error class c: Appendix B, Step 6(a)
			idx <- lc.list[[i]]  # index of locations in error class i
			y1.tmp <- y1.lc[[i]]
			y0.tmp <- y0.lc[[i]]
			s.y0 <- s[y0.tmp,]
			s.y1 <- s[y1.tmp,]
			mu.y0 <- mu.t[y0.tmp,]
			mu.y1 <- mu.t[y1.tmp,]

			###
			### Update sigma_c: Appendix B, Step 6(b)
			###

			if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				mh.star.sigma <- sum(dt2(s.y1,mu.y1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma.star,Q.star))
				mh.0.sigma <- sum(dt2(s.y1,mu.y1,z=1,Sigma,Q,i))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma,Q,i))
				if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					sigma[i] <- sigma.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$sigma[i] <- keep$sigma[i]+1
					keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				}
			}

			###
			### Update a_c: Appendix B, Step 6(c)
			###
			
			if(a.star[i]>0 & a.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				mh.star.a <- sum(dt2(s.y1,mu.y1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma.star,Q.star))
				mh.0.a <- sum(dt2(s.y1,mu.y1,z=1,Sigma,Q,i))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma,Q,i))
				if(exp(mh.star.a-mh.0.a)>runif(1)){
					a[i] <- a.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$a[i] <- keep$a[i]+1
					keep.tmp$a[i] <- keep.tmp$a[i]+1
				}
			}

			###
			### Update rho_c: Appendix B, Step 6(d)
			###
			
			if(rho.star[i]>0 & rho.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				mh.star.rho <- sum(dt2(s.y1,mu.y1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma.star,Q.star))
				mh.0.rho <- sum(dt2(s.y1,mu.y1,z=1,Sigma,Q,i))+
					sum(dt2(s.y0,mu.y0,z=0,Sigma,Q,i))
				if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					rho[i] <- rho.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$rho[i] <- keep$rho[i]+1
					keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				}
			}
		}

	
		#--------------------------------------------------------------------------
		# Update temporal process model parameters: Appendix B, Step 7
		#--------------------------------------------------------------------------

	 	t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		###	
		### Update fixed effects beta: Appendix B, Step 7(a)
		###
		
		b <- crossprod(X,(v-W%*%alpha))
		beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		###
		### Update basis coefficients alpha: Appendix B, Step 7(b)
		###
		
		b <- crossprod(W,(v-X%*%beta))
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		###
		### Update auxilliary variables v(t): Appendix B, Step 7(c)
		###

		linpred <- X%*%beta+W%*%alpha  # update linear predictor 
	 	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# End time of aux. variable update
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  

		###
		### Update variance of basis coefficients sigma_alpha^2: Appendix B, Step 7(d)
		###

		r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.alpha)
		q.tmp <- qW/2+q.alpha
		sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		diag(Sigma.alpha) <- sigma2.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)


		#--------------------------------------------------------------------------
		# Save samples: Appendix B, Step 8
		#--------------------------------------------------------------------------
						
		sigma.save[k,] <- sigma
		a.save[k,] <- a
		rho.save[k,] <- rho
		mu.save[,k] <- S.tilde[mu.j[h],2]
		theta.save[k] <- theta    
		sigma.mu.save[k] <- sigma.mu
		sigma.alpha.save[k] <- sigma2.alpha
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		v.save[,k] <- v
		m.save[k] <- m
		m.save.tmp <- m.save.tmp+m

	}  # End of MCMC loop

  	t.mcmc.end <- Sys.time()
  	
  	tune$sigma.mu <- tune$sigma.mu*s.sd
  	tune$sigma <- tune$sigma*s.sd
	tune$mu <- tune$mu*s.sd
  	sigma.save <- sigma.save*s.sd
  	sigma.mu.save <- sigma.mu.save*s.sd
	sigma.alpha.save <- sqrt(sigma.alpha.save)
	
	
  	#####################################################################
	### Write output
	#####################################################################

	# Ending values
	end <- list()
	end$sigma <- sigma*s.sd
	end$a <- a
	end$rho <- rho
	end$theta <- theta
	end$sigma.mu <- sigma.mu*s.sd
	end$mu <- S.tilde[mu.j,2]
	end$pi <- pi
	end$alpha <- alpha
	end$beta <- beta
	end$sigma.alpha <- sqrt(sigma2.alpha)
  
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu <- keep$mu/sum(m.save)
	keep$sigma <- keep$sigma/n.mcmc
	keep$a <- keep$a/n.mcmc
	keep$rho <- keep$rho/n.mcmc

	# Print acceptance rates
	cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	cat("\na acceptance rate:",round(keep$a,2))
	cat("\nrho acceptance rate:",round(keep$rho,2)) 

	# Print clocks
	cat(paste("\n\nEnd time:",Sys.time()))
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	cat(paste("\nTime per MCMC iteration:",
		round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	list(mu=mu.save,theta=theta.save,m=m.save,sigma.mu=sigma.mu.save,
		v=v.save,beta=beta.save,alpha=alpha.save,sigma.alpha=sigma.alpha.save,
		sigma=sigma.save,a=a.save,rho=rho.save,
		keep=keep,tune=tune,priors,J=J,start=start,priors=priors,n.mcmc=n.mcmc,end=end)
}