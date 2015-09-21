haulouts.2.mcmc <- function(s,y,X,W,S.tilde,sigma.alpha,
	priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (04 SEP 2015)
 	### See haulouts.1.sim.R to simulate data according to this model specification,
 	### and haulout.dp.mixture.2.pdf for the model description, model statement, and
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
  
  	truncnormsamp <- function(mu,sig2,low,high,nsamp){
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}

	# get.det <- function(X){
		# X[1,1]*X[2,2]-X[1,2]*X[2,1]
	# }

	# get.inv <- function(X,det){
		# det <- X[1,1]*X[2,2]-X[1,2]*X[2,1]
		# inv <- (1/det)*matrix(c(X[2,2],-X[2,1],-X[1,2],X[1,1]),2)		
		# list(inv=inv,det=det)
	# }
	
	get.Sigma <- function(sigma2,n.lc,Mix){  # get var-cov matrix and associated quantities
		Sigma <- lapply(1:n.lc,function(x) sigma2[x]*Mix[[x]])  # variance-covariance matrix
		det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		list(Sigma=Sigma,P=P,det=det)  # return list of lists
	}

	dmvnorm2 <- function(x,y,lc,Sigma,const=TRUE,log=TRUE,K=K){
	 	#Calculate log density of mixture normal distribution (Eq. 1)
		# browser()	
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
				# i <- 1
				idx <- lc.list[[i]]
				d.tmp <- d[idx,]
				P.tmp <- P[[i]]
				out[idx] <- exp(-0.5*rowSums((d.tmp%*%P.tmp)*d.tmp))+			
					exp(-0.5*rowSums((d.tmp%*%(K%*%P.tmp%*%t(K)))*d.tmp))
			}
			if(log) out <- log(out)  # log of mixture kernel
			if(const){  # add normalizing constant to kernel
				if(log){
					b <- log(0.5)+log(1/(2*pi*unlist(lapply(Sigma$det,function(x) x^(1/2)))))
					out <- out+b[lc]								
				} 
				if(!log){
					b <- 0.5/(2*pi*unlist(lapply(Sigma$det,function(x) x^(1/2))))
					out <- out*b[lc]			
				} 
			}
		}
		if(nrow(x)==0) out <- 0
		out
	}

	get.mh.mu <- function(x,mu.0,mu.0.star,S,ht,s,lc,z,Sigma,Sigma.z0,K=K){
		# browser()
		mu <- S[mu.0[x],3:4]
		mu.star <- S[mu.0.star[x],3:4]
		idx.0 <- which(ht==mu.0[x]&z==0)
		idx.1 <- which(ht==mu.0[x]&z==1)
		# U.0 <- U[ht[x],]
		# U.star <- U[ht.star[x],]
		mh.star <- sum(
			dmvnorm2(s[idx.1,],mu.star,lc[idx.1],Sigma,const=FALSE,log=TRUE),
			dmvnorm2(s[idx.0,],mu.star,lc[idx.0],Sigma.z0,const=FALSE,log=TRUE)) 
		mh.0 <-	sum(
			dmvnorm2(s[idx.1,],mu,lc[idx.1],Sigma,const=FALSE,log=TRUE),
			dmvnorm2(s[idx.0,],mu,lc[idx.0],Sigma.z0,const=FALSE,log=TRUE)) 
		exp(mh.star-mh.0)>runif(1)
	}
	
	get.ht <- function(x,y,z,lc,Sigma,Sigma.z0,const=FALSE,K=K){
	 	#Calculate log density of mixture normal distribution (Eq. 1)
		# browser()	
		if(z==1) P <- Sigma$P[[lc]]  # precision matrix if z==1
		if(z==0) P <- Sigma.z0$P[[lc]]  # precision matrix if z==0
		d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		out <- exp(-0.5*rowSums((d%*%P)*d))+exp(-0.5*rowSums((d%*%(K%*%P%*%t(K)))*d))
		log(out)  # log of mixture kernel
	}
	
	dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	#Calculate log density of mixture t distribution (Eq. 1)
		# browser()	
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
				P.tmp <- P[[i]]
				out[idx] <- (1+1/nu*(rowSums((d.tmp%*%P.tmp)*d.tmp)))^-((nu+2)/2) +
					(1+1/nu*(rowSums((d.tmp%*%(K%*%P.tmp%*%t(K)))*d.tmp)))^-((nu+2)/2)
			}	
			b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))
			if(log){
				out <- log(0.5)+log(out)+log(b[lc])
					# +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			} 
			if(!log){
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
  
	get.mh.mu <- function(x,mu.0,mu.0.star,S,ht,s,lc,z,Sigma,Sigma.z0){
		# browser()
		mu <- S[mu.0[x],3:4]
		mu.star <- S[mu.0.star[x],3:4]
		idx.0 <- which(ht==mu.0[x]&z==0)
		idx.1 <- which(ht==mu.0[x]&z==1)
		# U.0 <- U[ht[x],]
		# U.star <- U[ht.star[x],]
		mh.star <- sum(
			dmvt2(s[idx.1,],mu.star,lc[idx.1],Sigma,log=TRUE),
			dmvt2(s[idx.0,],mu.star,lc[idx.0],Sigma.z0,log=TRUE)) 
		mh.0 <-	sum(
			dmvt2(s[idx.1,],mu,lc[idx.1],Sigma,log=TRUE),
			dmvt2(s[idx.0,],mu,lc[idx.0],Sigma.z0,log=TRUE)) 
		exp(mh.star-mh.0)>runif(1)
	}

	get.ht <- function(x,y,z,lc,Sigma,Sigma.z0,nu=100,K=matrix(c(-1,0,0,1),2)){
	 	#Calculate log density of mixture normal distribution (Eq. 1)
		# browser()	
		if(z==1) P <- Sigma$P[[lc]]  # precision matrix if z==1
		if(z==0) P <- Sigma.z0$P[[lc]]  # precision matrix if z==0
		d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# out <- exp(-0.5*rowSums((d%*%P)*d))+exp(-0.5*rowSums((d%*%(K%*%P%*%t(K)))*d))
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
	v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	X.cross <- t(X)%*%X
	W.cross <- t(W)%*%W
	idx <- which(values(S.tilde)==1)
	S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))

	###
	### Standardize parameters
	###
# browser()
	cat("\nStandardizing variables....")
	s.scale <- max(apply(s,2,function(x) max(x)-min(x)))/6
	s.center <- apply(s,2,function(x) max(x)+min(x))/2
	s <- (s-matrix(s.center,nrow=T,ncol=2,byrow=TRUE))/s.scale
	tune$sigma <- tune$sigma/s.scale
	tune$sigma.mu <- tune$sigma.mu/s.scale
	tune$mu.0 <- tune$mu.0/s.scale
	S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.center,nrow=nrow(S.tilde),
		ncol=2,byrow=TRUE))/s.scale
	

	###
	### Starting values and priors
	###
# browser() 
	cat("\nGetting starting values....")
	beta <- matrix(start$beta,qX)
	alpha <- matrix(0,qW)
	# alpha <- start$alpha
	theta <- start$theta
	sigma.mu <- start$sigma.mu/s.scale
	pie <- start$pie
	z <- start$z

	# Observation model parameters (mixture t-distribution)
	lc <- as.numeric(priors$lc)  # Argos location quality class
	n.lc <- length(unique(lc))
	sigma <- priors$sigma/s.scale
	sigma2 <- sigma^2
	priors$sigma.mu.l <- priors$sigma.mu.l/s.scale
	priors$sigma.mu.u <- priors$sigma.mu.u/s.scale
	a <- priors$a
	# nu <- priors$nu
	rho <- priors$rho
	# K <- matrix(c(-1,0,0,1),2)  # mixture transformation matrix
	Mix <- lapply(1:n.lc, function(x)  # mixture matrix
		matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))
	Sigma <- get.Sigma(sigma2,n.lc,Mix)
	
	# Process model parameters
	H <- priors$H			
	Sigma.z0 <- get.Sigma(sigma2+sigma.mu^2,n.lc,Mix)
  	mu.beta <- matrix(0,qX,1)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	mu.alpha <- matrix(0,qW,1)
	Sigma.alpha <- diag(qW)*sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)

	# For updates of beta and alpha
	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	A.inv.beta <- solve(X.cross+Sigma.beta.inv)
	# A.inv <- solve(t(X)%*%X+Sigma.beta.inv)
	# A.inv <- solve(t(W)%*%W+Sigma.alpha.inv)

	y1 <- which(y==1)+T
	y0 <- which(y==0)+T
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	linpred <- X%*%beta+W%*%alpha


	###
	### Set up Dirichlet process mixture variables
	###

	ht <- match(start$ht,S.tilde[,2])
	tab <- table(ht)  # tabulate cluster membership
	m <- length(tab)  # number of clusters
	ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	mu.0 <- as.numeric(names(tab))  # idx of occupied clusters

  	
  	###
	### Create receptacles for output
	###
  
  	cat("\nCreating receptacles for output....")
	ht.save <- matrix(0,T,n.mcmc)
	theta.save <- numeric(n.mcmc)  # concentration parameter
	m.save <- numeric(n.mcmc)
	z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	v.save <- matrix(0,T+n,n.mcmc)
	alpha.save <- matrix(0,n.mcmc,qW)
	beta.save <- matrix(0,n.mcmc,qX)
	sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	# D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	keep <- list(mu.0=0,sigma.mu=0)
    
	###
	### Begin MCMC loop
	###

	t.v.update <- 0
	t.mcmc.start <- Sys.time()

  	cat("\nEntering MCMC Loop....\n")
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
		mh <- sapply(dup.idx,function(x)  # accepted proposals
			get.mh.mu(x,mu.0,mu.0.star,S.tilde,ht,s,lc,z,Sigma,Sigma.z0)) 
		keep$mu.0 <- keep$mu.0+sum(mh)
		keep.idx <- dup.idx[mh]
		mu.0[keep.idx] <- mu.0.star[keep.idx]

		# Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|S.tilde]
		idx <- sample(S.tilde[-mu.0,1],H-m,replace=FALSE)  # idx of new mu.0
		# sum(idx%in%ht)
	    
	    # Sample cluster assignment indicator, ht (Note: sampling order matters here)
		samp <- c(mu.0[ord],idx)  # sampling in order of decreasing membership
		# sum(duplicated(samp))
		ht <- sapply(1:T,function(x) sample(samp,1,prob= 
			exp(log(pie)+get.ht(s[x,],S.tilde[samp,3:4],	z[x],lc[x],Sigma,Sigma.z0))))
					
		# Tabulate cluster membership with base functions
		tab <- table(ht)  # tabulate cluster membership
		m <- length(tab)  # number of clusters
		ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
		mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		
	    # Stick-breaking process (Note: sampling order matters here)
		
		pad <- ifelse(m<H,H-m-1,0)
		tab.tmp <- c(tab[ord],rep(0,pad))  # membership in decreasing order
		eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # Sample theta (concentration parameter); See Gelman section 23.3
       	theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta
# theta <- 1
		 # theta.star <- rnorm(1,theta,0.25)
		    # if(theta.star>0 & theta.star<10){
		    	# mh.star.sigma <- sum(dbeta(eta[-H],1,theta.star,log=TRUE))
		    	# mh.0.sigma <- sum(dbeta(eta[-H],1,theta,log=TRUE))	    	
			    # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
		    	    # theta <- theta.star
					# # sigma.z0 <- sigma.z0.star
		        	# # Sigma.inv <- solve(sigma^2*diag(2))
			        # # keep$sigma <- keep$sigma+1
		    	# } 
		    # }

	
		###
		### Updates pertaining to wet/dry status of y and z
		###
# browser()
	 	t.v.start <- Sys.time()
	
		# Sample v(t.tilde) (auxilliary variable for y) 
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# Sample v(t) (auxilliary variable for s) 
		z1 <- which(z==1)
		z0 <- which(z==0)
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))

		#  Sample alpha (coefficients on basis expansion W)
	# browser()
		# b <- t(W)%*%(v-X%*%beta)
		b <- crossprod(W,(v-X%*%beta))
		# alpha <- A.inv.alpha%*%b+t(chol(A.inv.alpha))%*%matrix(rnorm(qW),qW,1)
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		#  Sample beta (coefficients on covariates influencing haul-out probability)
	  	# b <- t(X)%*%(v-W%*%alpha)  # +mu.beta%*%Sigma.beta.inv
		b <- crossprod(X,(v-W%*%alpha))
	  	# beta <- A.inv.beta%*%b+t(chol(A.inv.beta))%*%matrix(rnorm(qX),qX,1)
	  	beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")
	  	
	    # Sample z(t) (prediction of haul-out indicator variable for times t)
		linpred <- X%*%beta+W%*%alpha
		p <- pnorm(linpred[1:T,])
		# p.tmp1 <- p*dmvnorm2(s,S.tilde[ht,3:4],lc,Sigma,const=TRUE,log=FALSE,K=K)
		# p.tmp2 <- (1-p)*dmvnorm2(s,S.tilde[ht,3:4],lc,Sigma.z0,const=TRUE,log=FALSE,K=K)
		p.tmp1 <- p*dmvt2(s,S.tilde[ht,3:4],lc,Sigma,log=FALSE)
		p.tmp2 <- (1-p)*dmvt2(s,S.tilde[ht,3:4],lc,Sigma.z0,log=FALSE)
		p.tmp <- exp(log(p.tmp1)-log(p.tmp1+p.tmp2))	
# hist(p.tmp)
		z <- rbinom(T,1,p.tmp)
# z <- start$z

		
	    ###
	    ### Sample sigma.mu (disperson around homerange center)
	    ###
# browser()

	    sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			idx <- which(z==0)
			mu.0.tmp <- S.tilde[ht[idx],3:4]
			lc.tmp <- lc[idx]
			Sigma.z0.star <- get.Sigma(sigma2+sigma.mu.star^2,n.lc,Mix)
		    # mh.star.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0.star))
		    # mh.0.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0))
		    mh.star.sigma.mu <- sum(dmvt2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0.star,log=TRUE))
		    mh.0.sigma.mu <- sum(dmvt2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0,log=TRUE))
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
				Sigma.z0 <- Sigma.z0.star
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
		sigma.mu.save[k] <- sigma.mu
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		v.save[,k] <- v
		z.save[,k] <- z
		m.save[k] <- m
	}
  	
  	sigma.mu.save <- sigma.mu.save*s.scale
  	t.mcmc.end <- Sys.time()

	###
	### Write output
	###
	  
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu.0 <- keep$mu.0/sum(m.save)
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	cat(paste("\n\nEnd time:",Sys.time()))
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	cat(paste("\nTime per MCMC iteration:",
		round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	list(ht=ht.save,beta=beta.save,alpha=alpha.save,
		theta=theta.save,sigma.mu=sigma.mu.save,z=z.save,v=v.save,
	  	m=m.save,keep=keep,n.mcmc=n.mcmc)
}



# haulouts.2.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S.tilde,sigma.alpha,
	# priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# ###
 	# ### Brian M. Brost (04 SEP 2015)
 	# ### See haulouts.1.sim.R to simulate data according to this model specification,
 	# ### and haulout.dp.mixture.2.pdf for the model description, model statement, and
 	# ### full conditional distributions
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s; X.tilde=design matrix containing covariates
 	# ### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	# ### W.tilde=basis expansion for y; S.tilde=support Dirichlet process mixture 
 	# ### (i.e., haul-out sites); sigma.alpha=standard deviation of parameter 
 	# ### model for 'random' effects...
 	# ###
  
	# t.start <- Sys.time()

	# ###
	# ### Libraries and Subroutines
	# ###
  
	# # library(MCMCpack)  # for Dirichlet distribution functions
	# library(data.table)  # for tabulating and summing
	# # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# library(doParallel)  # for parallel processing
	# library(foreach)  # for parallel processing  
	# # library(mvtnorm)  # for multivariate normal density
	# # library(msm)  # for truncated normal density
  
  	# truncnormsamp <- function(mu,sig2,low,high,nsamp){
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }

	# # get.det <- function(X){
		# # X[1,1]*X[2,2]-X[1,2]*X[2,1]
	# # }

	# # get.inv <- function(X,det){
		# # det <- X[1,1]*X[2,2]-X[1,2]*X[2,1]
		# # inv <- (1/det)*matrix(c(X[2,2],-X[2,1],-X[1,2],X[1,1]),2)		
		# # list(inv=inv,det=det)
	# # }
	
	# get.Sigma <- function(sigma2,n.lc,Mix){  # get var-cov matrix and associated quantities
		# Sigma <- lapply(1:n.lc,function(x) sigma2[x]*Mix[[x]])  # variance-covariance matrix
		# det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		# P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			# matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		# list(Sigma=Sigma,P=P,det=det)  # return list of lists
	# }

	# get.mh.mu <- function(x,mu.0,mu.0.star,S,ht,s,lc,z,Sigma,Sigma.z0,K=K){
		# # browser()
		# mu <- S[mu.0[x],3:4]
		# mu.star <- S[mu.0.star[x],3:4]
		# idx.0 <- which(ht==mu.0[x]&z==0)
		# idx.1 <- which(ht==mu.0[x]&z==1)
		# # U.0 <- U[ht[x],]
		# # U.star <- U[ht.star[x],]
		# mh.star <- sum(
			# dmvnorm2(s[idx.1,],mu.star,lc[idx.1],Sigma,const=FALSE,log=TRUE),
			# dmvnorm2(s[idx.0,],mu.star,lc[idx.0],Sigma.z0,const=FALSE,log=TRUE)) 
		# mh.0 <-	sum(
			# dmvnorm2(s[idx.1,],mu,lc[idx.1],Sigma,const=FALSE,log=TRUE),
			# dmvnorm2(s[idx.0,],mu,lc[idx.0],Sigma.z0,const=FALSE,log=TRUE)) 
		# exp(mh.star-mh.0)>runif(1)
	# }
		
	# dmvnorm2 <- function(x,y,lc,Sigma,const=TRUE,log=TRUE,K=K){
	 	# #Calculate log density of mixture normal distribution (Eq. 1)
		# # browser()	
		# x <- matrix(x,,2,byrow=FALSE)
		# if(nrow(x)!=0){
			# if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# lc.idx <- sort(unique(lc))
			# n <- length(lc.idx)
			# lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# P <- Sigma$P[lc.idx]  # precision matrix
			# d <- x-y  
			# out <- numeric(nrow(d))
			# for(i in 1:n){  # calculate kernel for each Sigma
				# # i <- 1
				# idx <- lc.list[[i]]
				# d.tmp <- d[idx,]
				# P.tmp <- P[[i]]
				# out[idx] <- exp(-0.5*rowSums((d.tmp%*%P.tmp)*d.tmp))+			
					# exp(-0.5*rowSums((d.tmp%*%(K%*%P.tmp%*%t(K)))*d.tmp))
			# }
			# if(log) out <- log(out)  # log of mixture kernel
			# if(const){  # add normalizing constant to kernel
				# if(log){
					# b <- log(0.5)+log(1/(2*pi*unlist(lapply(Sigma$det,function(x) x^(1/2)))))
					# out <- out+b[lc]								
				# } 
				# if(!log){
					# b <- 0.5/(2*pi*unlist(lapply(Sigma$det,function(x) x^(1/2))))
					# out <- out*b[lc]			
				# } 
			# }
		# }
		# if(nrow(x)==0) out <- 0
		# out
	# }

	# get.ht <- function(x,y,z,lc,Sigma,Sigma.z0,const=FALSE,K=K){
	 	# #Calculate log density of mixture normal distribution (Eq. 1)
		# # browser()	
		# if(z==1) P <- Sigma$P[[lc]]  # precision matrix if z==1
		# if(z==0) P <- Sigma.z0$P[[lc]]  # precision matrix if z==0
		# d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# out <- exp(-0.5*rowSums((d%*%P)*d))+exp(-0.5*rowSums((d%*%(K%*%P%*%t(K)))*d))
		# log(out)  # log of mixture kernel
	# }

	# dmvt2 <- function(x,y,lc,Sigma,nu=100,K=K,const=TRUE,log=TRUE){
	 	# #Calculate log density of mixture t distribution (Eq. 1)
		# # browser()	
		# x <- matrix(x,,2,byrow=FALSE)
		# if(nrow(x)!=0){
			# if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# lc.idx <- sort(unique(lc))
			# n <- length(lc.idx)
			# lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# P <- Sigma$P[lc.idx]  # precision matrix
			# d <- x-y  
			# out <- numeric(nrow(d))
			# for(i in 1:n){  # calculate kernel for each Sigma
				# idx <- lc.list[[i]]
				# d.tmp <- d[idx,]
				# P.tmp <- P[[i]]
				# out[idx] <- (1+1/nu*(rowSums((d.tmp%*%P.tmp)*d.tmp)))^-((nu+2)/2) +
					# (1+1/nu*(rowSums((d.tmp%*%(K%*%P.tmp%*%t(K)))*d.tmp)))^-((nu+2)/2)
			# }	
			# if(log) out <- log(out)  # log of mixture kernel
			# if(const){  # add normalizing constant to kernel
				# b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))
				# if(log)	out <- log(0.5)+out+log(b[lc])
				# if(!log) out <- 0.5*out*b[lc]			
			# }
		# }
		# if(nrow(x)==0) out <- 0
		# out
	# }

  
	# ###
	# ###  Setup Variables 
	# ###
  
# # browser() 
	# T <- nrow(s)  # number of telemetry locations
	# n <- length(y)  # number of wet/dry observations
	# qX <- ncol(X)
	# qW <- ncol(W)
	# v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	# X.comb <- rbind(X,X.tilde)  # combined design matrix for updates on alpha, beta
	# W.comb <- rbind(W,W.tilde)  # combined design matrix for updates on alpha, beta 
	# idx <- which(values(S.tilde)==1)
	# S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))

	# ###
	# ### Standardize parameters
	# ###
# # browser()
	# cat("\nStandardizing variables....")
	# s.scale <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.center <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.center,nrow=T,ncol=2,byrow=TRUE))/s.scale
	# tune$sigma <- tune$sigma/s.scale
	# tune$sigma.mu <- tune$sigma.mu/s.scale
	# tune$mu.0 <- tune$mu.0/s.scale
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.center,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.scale
	

	# ###
	# ### Starting values and priors
	# ###
# # browser() 
	# cat("\nGetting starting values....")
	# beta <- matrix(start$beta,qX)
	# alpha <- matrix(0,qW)
	# # alpha <- start$alpha
	# theta <- start$theta
	# sigma.mu <- start$sigma.mu/s.scale
	# pie <- start$pie
	# z <- start$z

	# # Observation model parameters (mixture t-distribution)
	# lc <- as.numeric(priors$lc)  # Argos location quality class
	# n.lc <- length(unique(lc))
	# sigma <- priors$sigma/s.scale
	# sigma2 <- sigma^2
	# priors$sigma.mu.l <- priors$sigma.mu.l/s.scale
	# priors$sigma.mu.u <- priors$sigma.mu.u/s.scale
	# a <- priors$a
	# # nu <- priors$nu
	# rho <- priors$rho
	# K <- matrix(c(-1,0,0,1),2)  # mixture transformation matrix
	# Mix <- lapply(1:n.lc, function(x)  # mixture matrix
		# matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))
	# Sigma <- get.Sigma(sigma2,n.lc,Mix)
	
	# # Process model parameters
	# H <- priors$H			
	# Sigma.z0 <- get.Sigma(sigma2+sigma.mu^2,n.lc,Mix)
  	# mu.beta <- matrix(0,qX,1)
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# mu.alpha <- matrix(0,qW,1)
	# Sigma.alpha <- diag(qW)*sigma.alpha^2
  	# Sigma.alpha.inv <- solve(Sigma.alpha)

	# y1 <- which(y==1)+T
	# y0 <- which(y==0)+T
	# y1.sum <- length(y1)
	# y0.sum <- length(y0)
	# linpred <- X.comb%*%beta+W.comb%*%alpha


	# ###
	# ### Set up Dirichlet process mixture variables
	# ###

	# ht <- match(start$ht,S.tilde[,2])
	# tab <- table(ht)  # tabulate cluster membership
	# m <- length(tab)  # number of clusters
	# ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	# mu.0 <- as.numeric(names(tab))  # idx of occupied clusters

  	
  	# ###
	# ### Create receptacles for output
	# ###
  
  	# cat("\nCreating receptacles for output....")
	# ht.save <- matrix(0,T,n.mcmc)
	# theta.save <- numeric(n.mcmc)  # concentration parameter
	# m.save <- numeric(n.mcmc)
	# z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	# v.save <- matrix(0,T+n,n.mcmc)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# beta.save <- matrix(0,n.mcmc,qX)
	# sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	# # D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	# keep <- list(mu.0=0,sigma.mu=0)
    
	# ###
	# ### Begin MCMC loop
	# ###

  	# cat("\nEntering MCMC Loop....\n")
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# } 

		# ###
		# ### Dirichlet process parameters
		# ### Note: sampling order matters here. Cluster parameters must be 
		# ### sampled in same order as pie, i.e., sorted by decreasing membership
		# ###
	
		# # Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# # and Gelman et al. (2014), Section 23.3

		# # Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		   
		# # Note: sampling order does not matter here
# browser()			
		
		# # plot(S.tilde[,3:4],pch=19,cex=0.5,col="grey85")
		# # points(S.tilde[mu.0,3:4],pch=19,col=1)
		# mu.0.star <- sapply(mu.0,function(x)  # proposals for mu.0	
			# sample(S.tilde[,1],1,prob=dnorm(S.tilde[x,3],S.tilde[,3],tune$mu.0)*
			# dnorm(S.tilde[x,4],S.tilde[,4],tune$mu.0)))  
		# # points(S.tilde[mu.0.star,3:4],col=rgb(1,0,0,0.25),cex=0.5,pch=19)
		# dup.idx <- which(!duplicated(mu.0.star))  # exclude duplicate proposals
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu.0,mu.0.star,S.tilde,ht,s,lc,z,Sigma,Sigma.z0)) 
		# keep$mu.0 <- keep$mu.0+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu.0[keep.idx] <- mu.0.star[keep.idx]

		# # Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|S.tilde]
		# idx <- sample(S.tilde[-mu.0,1],H-m,replace=FALSE)  # idx of new mu.0
		# # sum(idx%in%ht)
	    
	    # # Sample cluster assignment indicator, ht (Note: sampling order matters here)
		# samp <- c(mu.0[ord],idx)  # sampling in order of decreasing membership
		# # sum(duplicated(samp))

	# # browser()			
		# ht <- sapply(1:T,function(x) sample(samp,1,prob= 
			# exp(log(pie)+get.ht(s[x,],S.tilde[samp,3:4],	z[x],lc[x],Sigma,Sigma.z0))))
					
		# # Tabulate cluster membership with base functions
		# tab <- table(ht)  # tabulate cluster membership
		# m <- length(tab)  # number of clusters
		# ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
		# mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		
	    # # Stick-breaking process (Note: sampling order matters here)
		# tab.tmp <- c(tab[ord],rep(0,H-m-1))  # membership in decreasing order
		# eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    # pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # # Sample theta (concentration parameter); See Gelman section 23.3
       	# theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta

		 # # theta.star <- rnorm(1,theta,0.25)
		    # # if(theta.star>0 & theta.star<10){
		    	# # mh.star.sigma <- sum(dbeta(eta[-H],1,theta.star,log=TRUE))
		    	# # mh.0.sigma <- sum(dbeta(eta[-H],1,theta,log=TRUE))	    	
			    # # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
		    	    # # theta <- theta.star
					# # # sigma.z0 <- sigma.z0.star
		        	# # # Sigma.inv <- solve(sigma^2*diag(2))
			        # # # keep$sigma <- keep$sigma+1
		    	# # } 
		    # # }

	
		# ###
		# ### Updates pertaining to wet/dry status of y and z
		# ###
# # browser()
		# # Sample v(t.tilde) (auxilliary variable for y) 
	  	# v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# # Sample v(t) (auxilliary variable for s) 
		# z1 <- which(z==1)
		# z0 <- which(z==0)
	  	# v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))

		# #  Sample alpha (coefficients on basis expansion W)
		# A.inv <- solve(t(W.comb)%*%W.comb+Sigma.alpha.inv)
		# b <- t(W.comb)%*%(v-X.comb%*%beta)
		# alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
		
		# #  Sample beta (coefficients on covariates influencing haul-out probability)
		# A.inv <- solve(t(X.comb)%*%X.comb+Sigma.beta.inv)
	  	# b <- t(X.comb)%*%(v-W.comb%*%alpha)  # +mu.beta%*%Sigma.beta.inv
	  	# beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	    # # Sample z(t) (prediction of haul-out indicator variable for times t)
		# linpred <- X.comb%*%beta+W.comb%*%alpha
		# p <- pnorm(linpred[1:T,])
		# p.tmp1 <- p*dmvnorm2(s,S.tilde[ht,3:4],lc,Sigma,const=TRUE,log=FALSE,K=K)
		# p.tmp2 <- (1-p)*dmvnorm2(s,S.tilde[ht,3:4],lc,Sigma.z0,const=TRUE,log=FALSE,K=K)
		# p.tmp <- exp(log(p.tmp1)-log(p.tmp1+p.tmp2))	
		# z <- rbinom(T,1,p.tmp)
# # z <- start$z

		 
	    # ###
	    # ### Sample sigma.mu (disperson around homerange center)
	    # ###
# # browser()

	    # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# idx <- which(z==0)
			# mu.0.tmp <- S.tilde[ht[idx],3:4]
			# lc.tmp <- lc[idx]
			# Sigma.z0.star <- get.Sigma(sigma2+sigma.mu.star^2,n.lc,Mix)
		    # mh.star.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0.star))
		    # mh.0.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,Sigma.z0))
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Sigma.z0 <- Sigma.z0.star
				# keep$sigma.mu <- keep$sigma.mu+1
	    	# } 
	    # }
# # sigma.mu <- start$sigma.mu
	
		# ###
		# ###  Save samples 
		# ###
		
# # browser()	
		# ht.save[,k] <- S.tilde[ht,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# v.save[,k] <- v
		# z.save[,k] <- z
		# m.save[k] <- m
	# }
  	
  	# sigma.mu.save <- sigma.mu.save*s.scale
  
	# ###
	# ### Write output
	# ###
	  
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu.0 <- keep$mu.0/sum(m.save)
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# list(ht=ht.save,beta=beta.save,alpha=alpha.save,
		# theta=theta.save,sigma.mu=sigma.mu.save,z=z.save,v=v.save,
	  	# m=m.save,keep=keep,n.mcmc=n.mcmc)
# }


# This version below works fine, but does not use the list form of Sigma (i.e., list
# containing the variance-covariance matrix, the precision matrix, and the determinant)

# haulouts.2.mcmc <- function(s,y,X,X.tilde,W,W.tilde,S.tilde,sigma.alpha,
	# priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# ###
 	# ### Brian M. Brost (04 SEP 2015)
 	# ### See haulouts.1.sim.R to simulate data according to this model specification,
 	# ### and haulout.dp.mixture.2.pdf for the model description, model statement, and
 	# ### full conditional distributions
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s; X.tilde=design matrix containing covariates
 	# ### influencing wet/dry status of ancillary data y; W=basis expansion for s; 
 	# ### W.tilde=basis expansion for y; S.tilde=support Dirichlet process mixture 
 	# ### (i.e., haul-out sites); sigma.alpha=standard deviation of parameter 
 	# ### model for 'random' effects...
 	# ###
  
	# t.start <- Sys.time()

	# ###
	# ### Libraries and Subroutines
	# ###
  
	# # library(MCMCpack)  # for Dirichlet distribution functions
	# library(data.table)  # for tabulating and summing
	# # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# library(doParallel)  # for parallel processing
	# library(foreach)  # for parallel processing  
	# # library(mvtnorm)  # for multivariate normal density
	# # library(msm)  # for truncated normal density
  
  	# truncnormsamp <- function(mu,sig2,low,high,nsamp){
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }

	# get.det <- function(X){  
		# X[1,1]*X[2,2]-X[1,2]*X[2,1]
	# }

	# get.inv <- function(X,det){
		# det <- X[1,1]*X[2,2]-X[1,2]*X[2,1]
		# (1/det)*matrix(c(X[2,2],-X[2,1],-X[1,2],X[1,1]),2)		
	# }

# get.mh.mu <- function(x,mu.0,mu.0.star,S,ht,s,lc,z,Sig.i,Sig.d,Sig.z0.i,Sig.z0.d){
		# # browser()
		# # x <- 2
		# mu <- S[mu.0[x],3:4]
		# mu.star <- S[mu.0.star[x],3:4]
		# idx.0 <- which(ht==mu.0[x]&z==0)
		# idx.1 <- which(ht==mu.0[x]&z==1)
		# # U.0 <- U[ht[x],]
		# # U.star <- U[ht.star[x],]
		# mh.star <- sum(
			# dmvnorm2(s[idx.1,],mu.star,lc[idx.1],Sig.i,Sig.d,const=TRUE,log=TRUE),
			# dmvnorm2(s[idx.0,],mu.star,lc[idx.0],Sig.z0.i,Sig.z0.d,const=TRUE,log=TRUE)) 
		# mh.0 <-	sum(
			# dmvnorm2(s[idx.1,],mu,lc[idx.1],Sig.i,Sig.d,const=TRUE,log=TRUE),
			# dmvnorm2(s[idx.0,],mu,lc[idx.0],Sig.z0.i,Sig.z0.d,const=TRUE,log=TRUE)) 
		# exp(mh.star-mh.0)>runif(1)
	# }
		
	# dmvnorm2 <- function(x,y,lc,Sig.i,Sig.d,const=TRUE,log=TRUE){
	 	# #Calculate log density of mixture normal distribution (Eq. 1)
		# # browser()	
		# x <- matrix(x,,2,byrow=FALSE)
		# if(nrow(x)!=0){
			# if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# K <- matrix(c(-1,0,0,1),2)  # mixture transformation matrix
			# lc.idx <- sort(unique(lc))
			# n <- length(lc.idx)
			# lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# P <- Sig.i[lc.idx]  # precision matrix
			# d <- x-y  
			# out <- numeric(nrow(d))
			# for(i in 1:n){  # calculate kernel for each Sigma
				# # i <- 1
				# idx <- lc.list[[i]]
				# d.tmp <- d[idx,]
				# P.tmp <- P[[i]]
				# out[idx] <- exp(-0.5*rowSums((d.tmp%*%P.tmp)*d.tmp))+			
					# exp(-0.5*rowSums((d.tmp%*%(K%*%P.tmp%*%t(K)))*d.tmp))
			# }
			# if(log) out <- log(out)  # log of mixture kernel
			# if(const){  # add normalizing constant to kernel
				# if(log){
					# b <- log(0.5)+log(1/(2*pi*unlist(lapply(Sig.d,function(x) x^(1/2)))))
					# out <- out+b[lc]								
				# } 
				# if(!log){
					# b <- 0.5/(2*pi*unlist(lapply(Sig.d,function(x) x^(1/2))))
					# out <- out*b[lc]			
				# } 
			# }
		# }
		# if(nrow(x)==0) out <- 0
		# out
	# }


  
	# ###
	# ###  Setup Variables 
	# ###
  
# # browser() 
	# T <- nrow(s)  # number of telemetry locations
	# n <- length(y)  # number of wet/dry observations
	# qX <- ncol(X)
	# qW <- ncol(W)
	# v <- numeric(n+T)  # auxilliary variable for continuous haul-out process
	# X.comb <- rbind(X,X.tilde)  # combined design matrix for updates on alpha, beta
	# W.comb <- rbind(W,W.tilde)  # combined design matrix for updates on alpha, beta 
	# idx <- which(values(S.tilde)==1)
	# S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))

	# ###
	# ### Standardize parameters
	# ###
# # browser()
	# cat("\nStandardizing variables....")
	# s.scale <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.center <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.center,nrow=T,ncol=2,byrow=TRUE))/s.scale
	# tune$sigma <- tune$sigma/s.scale
	# tune$sigma.mu <- tune$sigma.mu/s.scale
	# tune$mu.0 <- tune$mu.0/s.scale
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.center,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.scale
	

	# ###
	# ### Starting values and priors
	# ###
# # browser() 
	# beta <- matrix(start$beta,qX)
	# alpha <- matrix(0,qW)
	# # alpha <- start$alpha
	# theta <- start$theta
	# sigma.mu <- start$sigma.mu/s.scale
	# pie <- start$pie
	# z <- start$z

	# # Observation model parameters (mixture t-distribution)
	# lc <- as.numeric(priors$lc)  # Argos location quality class
	# n.lc <- length(unique(lc))
	# sigma <- priors$sigma/s.scale
	# priors$sigma.mu.l <- priors$sigma.mu.l/s.scale
	# priors$sigma.mu.u <- priors$sigma.mu.u/s.scale
	# a <- priors$a
	# # nu <- priors$nu
	# rho <- priors$rho
			
	# H <- priors$H

# # browser()	
	# Mix <- lapply(1:n.lc, function(x)  # mixture matrix
		# matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))

	# Sigma <- lapply(1:n.lc, function(x) sigma[x]^2*Mix[[x]])
	# Sigma.inv <- lapply(Sigma,solve)
	# Sigma.det <- lapply(Sigma,det)

  	# mu.beta <- matrix(0,qX,1)
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# mu.alpha <- matrix(0,qW,1)
	# Sigma.alpha <- diag(qW)*sigma.alpha^2
  	# Sigma.alpha.inv <- solve(Sigma.alpha)

	# # sigma.z0 <- sqrt(sigma^2+sigma.mu^2)
	# Sigma.z0 <- lapply(1:n.lc,function(x) Sigma[[x]]+sigma.mu^2*Mix[[x]])	
	# Sigma.z0.inv <- lapply(Sigma.z0,solve)
	# Sigma.z0.det <- lapply(Sigma.z0,det)

	# y1 <- which(y==1)+T
	# y0 <- which(y==0)+T
	# y1.sum <- length(y1)
	# y0.sum <- length(y0)
	# linpred <- X.comb%*%beta+W.comb%*%alpha

# # browser()	
	# ###
	# ### Set up Dirichlet process mixture variables
	# ###

	# ht <- match(start$ht,S.tilde[,2])
	# tab <- table(ht)  # tabulate cluster membership
	# m <- length(tab)  # number of clusters
	# ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
	# mu.0 <- as.numeric(names(tab))  # idx of occupied clusters

  	
  	# ###
	# ### Create receptacles for output
	# ###
  
	# ht.save <- matrix(0,T,n.mcmc)
	# theta.save <- numeric(n.mcmc)  # concentration parameter
	# m.save <- numeric(n.mcmc)
	# z.save <- matrix(0,T,n.mcmc)  # haul-out indicator variable
	# v.save <- matrix(0,T+n,n.mcmc)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# beta.save <- matrix(0,n.mcmc,qX)
	# sigma.mu.save <- numeric(n.mcmc)  # dispersion about haul-out site
	# # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of parameter model
	# D.bar.save <- numeric(n.mcmc)  # D.bar for DIC calculation

	# keep <- list(mu.0=0,sigma=0,sigma.mu=0)
    
	# ###
	# ### Begin MCMC loop
	# ###
  
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# } 

		# ###
		# ### Dirichlet process parameters
		# ### Note: sampling order matters here. Cluster parameters must be 
		# ### sampled in same order as pie, i.e., sorted by decreasing membership
		# ###
	
		# # Update follows the blocked Gibbs sampler of Ishwaran and James (2001)
		# # and Gelman et al. (2014), Section 23.3

		# # Sample 'occupied' mu.0 (true location of clusters with non-zero membership)		   
		# # Note: sampling order does not matter here
# # browser()			
		# # plot(S.tilde[,3:4],pch=19,cex=0.5,col="grey85")
		# # points(S.tilde[mu.0,3:4],pch=19,col=1)
		# mu.0.star <- sapply(mu.0,function(x)  # proposals for mu.0	
			# sample(S.tilde[,1],1,prob=dnorm(S.tilde[x,3],S.tilde[,3],tune$mu.0)*
			# dnorm(S.tilde[x,4],S.tilde[,4],tune$mu.0)))  
		# # points(S.tilde[mu.0.star,3:4],col=rgb(1,0,0,0.25),cex=0.5,pch=19)
		# dup.idx <- which(!duplicated(mu.0.star))  # exclude duplicate proposals
		
	# # browser()	
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu.0,mu.0.star,S.tilde,ht,s,lc,z,Sigma.inv,Sigma.det,Sigma.z0.inv,
			# Sigma.z0.det)) 

		# keep$mu.0 <- keep$mu.0+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu.0[keep.idx] <- mu.0.star[keep.idx]

		# # Sample 'unoccupied' mu.0 (clusters with zero membership) from prior, [m.0|S.tilde]
		# idx <- sample(S.tilde[-mu.0,1],H-m,replace=FALSE)  # idx of new mu.0
		# # sum(idx%in%ht)
	    
	    # # Sample cluster assignment indicator, ht (Note: sampling order matters here)
		# samp <- c(mu.0[ord],idx)  # sampling in order of decreasing membership
		# # sum(duplicated(samp))

	# get.ht <- function(x,y,z,lc,Sig.i,Sig.d,Sig.z0.i,Sig.z0.d,const=FALSE){
	 	# #Calculate log density of mixture normal distribution (Eq. 1)
		# # browser()	
		# K <- matrix(c(-1,0,0,1),2)  # mixture transformation matrix
		# if(z==1){
			# Sig.i <- Sig.i[[lc]]  # precision matrix if z==1
			# Sig.d <- Sig.d[[lc]]  # corresponding determinant
		# }
		# if(z==0){
			# Sig.i <- Sig.z0.i[[lc]]  # precision matrix if z==0
			# Sig.d <- Sig.z0.d[[lc]]	 # corresponding determinant
		# }
		# d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# out <- exp(-0.5*rowSums((d%*%Sig.i)*d))+	exp(-0.5*rowSums((d%*%(K%*%Sig.i%*%t(K)))*d))
		# out <- log(out)  # log of mixture kernel
		# # if(const){  # add normalizing constant to kernel
			# # b <- log(0.5)+log(1/(2*pi*Sig.d^(1/2)))
			# # out <- out+b
		# # }
		# out
	# }

# # browser()			
	
		# ht <- sapply(1:T,function(x) sample(samp,1,prob= 
			# exp(log(pie)+get.ht(s[x,],S.tilde[samp,3:4],
			# z[x],lc[x],Sigma.inv,Sigma.det,Sigma.z0.inv,Sigma.z0.det))))
					
		# # Tabulate cluster membership with base functions
		# tab <- table(ht)  # tabulate cluster membership
		# m <- length(tab)  # number of clusters
		# ord <- order(tab,decreasing=TRUE)  # clusters ordered by membership
		# mu.0 <- as.numeric(names(tab))  # idx of occupied clusters
		
	    # # Stick-breaking process (Note: sampling order matters here)
		# tab.tmp <- c(tab[ord],rep(0,H-m-1))  # membership in decreasing order
		# eta <- c(rbeta(H-1,1+tab.tmp,theta+T-cumsum(tab.tmp)),1)  # stick-breaking weights
	    # pie <- eta*c(1,cumprod((1-eta[-H])))  # mixture component probabilities

	    # # Sample theta (concentration parameter); See Gelman section 23.3
       	# theta <- rgamma(1,priors$r+H-1,priors$q-sum(log(1-eta[-H])))  
# theta <- start$theta

		 # # theta.star <- rnorm(1,theta,0.25)
		    # # if(theta.star>0 & theta.star<10){
		    	# # mh.star.sigma <- sum(dbeta(eta[-H],1,theta.star,log=TRUE))
		    	# # mh.0.sigma <- sum(dbeta(eta[-H],1,theta,log=TRUE))	    	
			    # # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
		    	    # # theta <- theta.star
					# # # sigma.z0 <- sigma.z0.star
		        	# # # Sigma.inv <- solve(sigma^2*diag(2))
			        # # # keep$sigma <- keep$sigma+1
		    	# # } 
		    # # }

	
		# ###
		# ### Updates pertaining to wet/dry status of y and z
		# ###
# # browser()
		# # Sample v(t.tilde) (auxilliary variable for y) 
	  	# v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# # Sample v(t) (auxilliary variable for s) 
		# z1 <- which(z==1)
		# z0 <- which(z==0)
	  	# v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))

		# #  Sample alpha (coefficients on basis expansion W)
		# A.inv <- solve(t(W.comb)%*%W.comb+Sigma.alpha.inv)
		# b <- t(W.comb)%*%(v-X.comb%*%beta)
		# alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
		
		# #  Sample beta (coefficients on covariates influencing haul-out probability)
		# A.inv <- solve(t(X.comb)%*%X.comb+Sigma.beta.inv)
	  	# b <- t(X.comb)%*%(v-W.comb%*%alpha)  # +mu.beta%*%Sigma.beta.inv
	  	# beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	    # # Sample z(t) (prediction of haul-out indicator variable for times t)
		# linpred <- X.comb%*%beta+W.comb%*%alpha
		# p <- pnorm(linpred[1:T,])
		# p.tmp1 <- p*dmvnorm2(s,S.tilde[ht,3:4],lc,Sigma.inv,Sigma.det,const=TRUE,log=FALSE)
		# p.tmp2 <- (1-p)*dmvnorm2(s,S.tilde[ht,3:4],lc,
			# Sigma.z0.inv,Sigma.z0.det,const=TRUE,log=FALSE)
		# p.tmp <- exp(log(p.tmp1)-log(p.tmp1+p.tmp2))	
		# z <- rbinom(T,1,p.tmp)
# # z <- start$z

		 
	    # ###
	    # ### Sample sigma.mu (disperson around homerange center)
	    # ###
# # browser()

	    # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# idx <- which(z==0)
			# mu.0.tmp <- S.tilde[ht[idx],3:4]
			# lc.tmp <- lc[idx]
			# Sigma.z0.star <- lapply(1:n.lc,function(x) Sigma[[x]]+sigma.mu.star^2*Mix[[x]])	
			# Sigma.z0.inv.star <- lapply(Sigma.z0.star,solve)
			# Sigma.z0.det.star <- lapply(Sigma.z0.star,det)
		    # mh.star.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,
		    	# Sigma.z0.inv.star,Sigma.z0.det.star,const=TRUE,log=TRUE))
		    # mh.0.sigma.mu <- sum(dmvnorm2(s[idx,],mu.0.tmp,lc.tmp,
		    	# Sigma.z0.inv,Sigma.z0.det,const=TRUE,log=TRUE))
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Sigma.z0.inv <- Sigma.z0.inv.star
				# Sigma.z0.det <- Sigma.z0.det.star
		        # keep$sigma.mu <- keep$sigma.mu+1
	    	# } 
	    # }
# # sigma.mu <- start$sigma.mu
	
		# ###
		# ###  Save samples 
		# ###
		
# # browser()	
		# ht.save[,k] <- S.tilde[ht,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# v.save[,k] <- v
		# z.save[,k] <- z
		# m.save[k] <- m
	# }
  	
  	# sigma.mu.save <- sigma.mu.save*s.scale
  
	# ###
	# ### Write output
	# ###
	  
	# keep$sigma <- keep$sigma/n.mcmc
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu.0 <- keep$mu.0/sum(m.save)
	# cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# cat(paste("\nmu.0 acceptance rate:",round(keep$mu.0,2))) 
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# list(ht=ht.save,beta=beta.save,alpha=alpha.save,
		# theta=theta.save,sigma.mu=sigma.mu.save,z=z.save,v=v.save,
	  	# m=m.save,keep=keep,n.mcmc=n.mcmc)
# }