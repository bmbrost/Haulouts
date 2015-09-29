haulouts.3.mcmc <- function(s,y,X,W,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (29 SEP 2015)
 	### See haulouts.3.sim.R to simulate data according to this model specification,
 	### and haulouts.3.pdf for the model description, model statement, and
 	### full conditional distributions
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s; X.tilde=design matrix containing covariates
 	### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	### W.tilde=basis expansion for y; 
 	
#	U
 
 	# S.tilde=support Dirichlet process mixture 
 	### (i.e., haul-out sites); sigma.alpha=standard deviation of parameter 
 	### model for 'random' effects...
 	###
  
	t.start <- Sys.time()
	cat(paste("Start time:",t.start,"\n"))

	###
	### Libraries and Subroutines
	###
  
	# library(MCMCpack)  # for Dirichlet distribution functions
	# library(data.table)  # for tabulating and summing
	# library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# library(doParallel)  # for parallel processing
	# library(foreach)  # for parallel processing  
	# library(mvtnorm)  # for multivariate normal density
	# library(msm)  # for truncated normal density
  
  	truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}

	get.Sigma <- function(sigma2,n.lc,Mix){  # get var-cov matrix, determinant, and inverse
		Sigma <- lapply(1:n.lc,function(x) sigma2[x]*Mix[[x]])  # variance-covariance matrix
		det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		list(Sigma=Sigma,P=P,det=det)  # return list of lists
	}

	dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# Calculate density of mixture t distribution
		x <- matrix(x,,2,byrow=FALSE)
		if(nrow(x)!=0){
			if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			lc.idx <- sort(unique(lc))
			n <- length(lc.idx)
			lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			P <- Sigma$P[lc.idx]  # precision matrix
			d <- x-y  
			out <- numeric(nrow(d))
			for(i in 1:n){  # calculate kernel for each Sigma
				idx <- lc.list[[i]]
				d.tmp <- d[idx,]
				out[idx] <- (1+1/nu*(rowSums((d.tmp%*%P[[i]])*d.tmp)))^-((nu+2)/2) +
					(1+1/nu*(rowSums((d.tmp%*%(K%*%P[[i]]%*%t(K)))*d.tmp)))^-((nu+2)/2)
			}	
			b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))  # determinant
			if(log){  # log density
				out <- log(0.5)+log(out)+log(b[lc])
					# +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			} 
			if(!log){  # density
				out <- 0.5*out*b[lc]
					# *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
			} 
		}
		if(nrow(x)==0) out <- 0
		out
	}

	# Test dmvt2 function
	# test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# test2 <- log(0.5)+
		# log(sapply(1:T,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# df=100,log=FALSE))+sapply(1:T,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# plot(test1,test2)
	# summary(test2-test1)
	# lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# test2 <- 0.5*
		# (sapply(1:T,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# df=100,log=FALSE))+sapply(1:T,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# plot(test1,test2)
	# summary(test2/test1)
	# gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	get.mh.mu <- function(x,mu.0,mu.0.star,S,ht,s,lc,z,Sigma,Sigma.z0,U,gamma){
		# Accept/reject proposals for mu.0
		mu <- S[mu.0[x],3:4]  # location of current clusters mu.0
		mu.star <- S[mu.0.star[x],3:4]  # location of proposal clusters mu.0.star
		idx.0 <- which(ht==mu.0[x]&z==0)  # obs. associated with mu.0 and z=0
		idx.1 <- which(ht==mu.0[x]&z==1)  # obs. associated with mu.0 and z=1
		U.0 <- U[mu.0[x],]  # resource selection covariates for mu.0
		U.star <- U[mu.0.star[x],]  # resource selection covariates for mu.0.star
		mh.star <- sum(  # numerator of Metropolis-Hastings ratio
			dmvt2(s[idx.1,],mu.star,lc[idx.1],Sigma,log=TRUE),
			dmvt2(s[idx.0,],mu.star,lc[idx.0],Sigma.z0,log=TRUE))+U.star%*%gamma
		mh.0 <-	sum(  # denominator of Metropolis-Hastings ratio
			dmvt2(s[idx.1,],mu,lc[idx.1],Sigma,log=TRUE),
			dmvt2(s[idx.0,],mu,lc[idx.0],Sigma.z0,log=TRUE))+U.0%*%gamma
		exp(mh.star-mh.0)>runif(1)  # Accept or reject
	}

	get.ht <- function(x,y,z,lc,Sigma,Sigma.z0,nu=100,K=matrix(c(-1,0,0,1),2)){
		# For sampling ht, assignment of obs. to current set of clusters mu.0
		if(z==1) P <- Sigma$P[[lc]]  # precision matrix when z=1 
		if(z==0) P <- Sigma.z0$P[[lc]]  # precision matrix when z=0
		d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		out <- (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2) +
			(1+1/nu*(rowSums((d%*%(K%*%P%*%t(K)))*d)))^-((nu+2)/2)
		log(out)  # log of mixture kernel
	}


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
	y1 <- which(y==1)+T
	y0 <- which(y==0)+T
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	X.cross <- t(X)%*%X  # cross product of X
	W.cross <- t(W)%*%W  # cross product of W
	idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# information in S.tilde; note that mu.0 below references row idx in S.tilde


	###
	### Standardize parameters
	###

	cat("\nStandardizing variables....")
# browser()	
	# Center and scale s and S.tilde
	s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	s <- (s-matrix(s.mean,nrow=T,ncol=2,byrow=TRUE))/s.sd
	S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		ncol=2,byrow=TRUE))/s.sd
	
	# Center and scale tuning parameters
	tune$sigma.mu <- tune$sigma.mu/s.sd
	tune$mu.0 <- tune$mu.0/s.sd
		

	###
	### Priors
	###
	
	# Observation model
	lc <- as.numeric(priors$lc)  # Argos location quality class
	n.lc <- length(unique(lc))  # number of location classes
	sigma <- priors$sigma/s.sd  # variance component
	sigma2 <- sigma^2  
	a <- priors$a  # variance component
	rho <- priors$rho  # variance component
	Mix <- lapply(1:n.lc, function(x)  # mixture observation model matrix
		matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))
	Sigma <- get.Sigma(sigma2,n.lc,Mix)  # Sigmas, determinants, and inverses
	Sigma.z0 <- get.Sigma(sigma2+sigma.mu^2,n.lc,Mix)  # Sigmas when at sea, i.e., z=0

	# Process model
	H <- priors$H  # maximum number of clusters per truncation approximation
	mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	A.inv.beta <- solve(X.cross+Sigma.beta.inv)
	mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	Sigma.gamma.inv <- solve(Sigma.gamma)

	# Center and scale priors
	priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu


	###
	### Starting values
	###

	cat("\nGetting starting values....")
	beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	linpred <- X%*%beta+W%*%alpha
	gamma <- matrix(start$gamma/U.sd,qU)  # haul-out location RSF coefficients
	gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	theta <- start$theta  # DP concentration parameter
	pie <- start$pie  # stick-breaking probability /pi
	z <- start$z  # haul-out status; z=1: dry, z=0: wet
	sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)
  	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)


	###
	### Set up Dirichlet process mixture variables
	###

	ht <- match(start$ht,S.tilde[,2])  # ht corresponds to row idx of S.tilde
	tab <- table(ht)  # tabulate cluster membership
	m <- length(tab)  # number of clusters
	ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		# mu.0 corresponds to row idx of S.tilde
  	

  	###
	### Create receptacles for output
	###
  
  	cat("\nCreating receptacles for output....")
	beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	ht.save <- matrix(0,T,n.mcmc)  # DP cluster assignment indicator
	theta.save <- numeric(n.mcmc)  # DP concentration parameter
	m.save <- numeric(n.mcmc)  # number of clusters
	v.save <- matrix(0,T+n,n.mcmc)  # auxiliary variable for temporal haul-out process
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

    
	###
	### Begin MCMC loop
	###

	keep <- list(mu.0=0,sigma.mu=0,gamma=0)  # number of proposals accepted for M-H updates
	t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	t.mcmc.start <- Sys.time()  # timing MCMC iterations
# browser()
  	cat("\nEntering MCMC Loop....\n")
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) {  # Monitor the appropriateness of H, the truncation approximation
	    	cat(k,"");flush.console()	
    		plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	} 

		###
		### Dirichlet process parameters
		###

		# Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# and Gelman et al. (2014), Section 23.3
		
		# Note: sampling order matters here. Cluster parameters must be 
		# sampled in same order as pie, i.e., sorted by decreasing membership
		
		# Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		   
		mu.0.star <- sapply(mu.0,function(x)  # proposals for mu.0	
			sample(S.tilde[,1],1,prob=dnorm(S.tilde[x,3],S.tilde[,3],tune$mu.0)*
			dnorm(S.tilde[x,4],S.tilde[,4],tune$mu.0)))  
		dup.idx <- which(!duplicated(mu.0.star))  # exclude duplicate proposals
		mh <- sapply(dup.idx,function(x)  # accepted proposals
			get.mh.mu(x,mu.0,mu.0.star,S.tilde,ht,s,lc,z,Sigma,Sigma.z0,U,gamma)) 
		keep$mu.0 <- keep$mu.0+sum(mh)
		keep.idx <- dup.idx[mh]
		mu.0[keep.idx] <- mu.0.star[keep.idx]  # update mu.0

		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior 
		p <- exp(U[-mu.0,]%*%gamma)
		idx <- sample(S.tilde[-mu.0,1],H-m,replace=FALSE,prob=p)  # idx of new mu.0

	    # Sample cluster assignment indicator
		samp <- c(mu.0[ord],idx)  # sampling in order of decreasing membership
		ht <- sapply(1:T,function(x) sample(samp,1,prob= 
			exp(log(pie)+get.ht(s[x,],S.tilde[samp,3:4],z[x],lc[x],Sigma,Sigma.z0))))
					
		# Tabulate cluster membership with base functions
		tab <- table(ht)  # tabulate cluster membership
		m <- length(tab)  # number of clusters
		ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
		mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		
	    # Stick-breaking process
		pad <- ifelse(m<H,H-m-1,0)
		tab.tmp <- c(tab[ord],rep(0,pad))  # membership in decreasing order
		eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # Sample theta (concentration parameter)
	    
	    # See Gelman section 23.3, Ishwaran and Zarepour (2000)
       	# theta <- rgamma(1,priors$r.theta+H-1,priors$q.theta-sum(log(1-eta[-H])))  

		# Escobar and West (1995) and West (1997?) white paper
		tmp <- rbeta(1,theta+1,T)
		c <- priors$r.theta
		d <- priors$q.theta
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
		### Updates pertaining to temporal haul-out process, i.e., wet/dry status of y and z
		###

	 	t.v.start <- Sys.time()  # timing update of auxiliary variable
	
		# Sample v(t.tilde) (auxilliary variable for wet/dry status y) 
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# Sample v(t) (auxilliary variable for telemetry observations s) 
		z1 <- which(z==1)
		z0 <- which(z==0)
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))

		#  Sample alpha (random effects for temporal haul-out process)
# browser()
		b <- crossprod(W,(v-X%*%beta))
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		#  Sample beta (temporal haul-out process coefficients; fixed effects)
		b <- crossprod(X,(v-W%*%alpha))
	  	beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")
	  	
	    # Sample z (haul-out indicator variable for telemetry locations at times t)
		linpred <- X%*%beta+W%*%alpha
		p <- pnorm(linpred[1:T,])
		p.tmp1 <- p*dmvt2(s,S.tilde[ht,3:4],lc,Sigma,log=FALSE)
		p.tmp2 <- (1-p)*dmvt2(s,S.tilde[ht,3:4],lc,Sigma.z0,log=FALSE)
		p.tmp <- exp(log(p.tmp1)-log(p.tmp1+p.tmp2))	
		z <- rbinom(T,1,p.tmp)


		###
	    ### Sample gamma (haul-out location RSF coefficients)
	    ###

		gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		gamma.int.star <- log(sum(exp(U%*%gamma.star)))
		mh.star.gamma <- sum(U[mu.0,]%*%gamma.star)-m*gamma.int.star+
			sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))
		mh.0.gamma <- sum(U[mu.0,]%*%gamma)-m*gamma.int+
			sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))
		if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    gamma <- gamma.star
	        gamma.int <- gamma.int.star
	        keep$gamma <- keep$gamma+1
    	} 

		
	    ###
	    ### Sample sigma.mu (homerange dispersion parameter)
	    ###

	    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			idx <- which(z==0)
			mu.0.tmp <- S.tilde[ht[idx],3:4]
			lc.tmp <- lc[idx]
			Sigma.z0.star <- get.Sigma(sigma2+sigma.mu.star^2,n.lc,Mix)
		    mh.star.sigma.mu <- sum(dmvt2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0.star,log=TRUE))
		    mh.0.sigma.mu <- sum(dmvt2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0,log=TRUE))
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
				Sigma.z0 <- Sigma.z0.star
				keep$sigma.mu <- keep$sigma.mu+1
	    	} 
	    }


		###
		###  Sample sigma.alpha (standard deviation of normal prior on alpha) 
		###

		r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		q.tmp <- qW/2+priors$q.sigma.alpha
		sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		diag(Sigma.alpha) <- sigma2.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)

			
		###
		###  Save samples 
		###
		
		ht.save[,k] <- S.tilde[ht,2]
		theta.save[k] <- theta    
		sigma.mu.save[k] <- sigma.mu
		sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		gamma.save[k,] <- gamma
		v.save[,k] <- v
		z.save[,k] <- z
		m.save[k] <- m
	}
  	
  	sigma.mu.save <- sigma.mu.save*s.sd
  	t.mcmc.end <- Sys.time()

	###
	### Write output
	###
	  
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu.0 <- keep$mu.0/sum(m.save)
	keep$gamma <- keep$gamma/n.mcmc
	cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat(paste("\n\nEnd time:",Sys.time()))
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	cat(paste("\nTime per MCMC iteration:",
		round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	list(beta=beta.save,alpha=alpha.save,gamma=gamma.save,
		ht=ht.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		sigma.mu=sigma.mu.save,sigma.alpha=sigma.alpha.save,
		keep=keep,n.mcmc=n.mcmc)
}
