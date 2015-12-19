haulouts.5.mcmc <- function(s,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (08 DEC 2015)
 	### Cleaned up version of haulouts.3.mcmc.R
 	###
 	
 	###
 	### Function arguments: s=telemetry locations, y=ancillary data source containing
 	### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s and wet/dry status of ancillary data y;  
 	### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	### the location of haul-out sites; S.tilde=support of haul-out sites 
 	###

	t.start <- Sys.time()
	cat(paste("Start time:",t.start,"\n"))

	#####################################################################
	### Libraries and Subroutines
	#####################################################################
  
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
	
	tkern <- function(d,P,nu=100){  # kernel of t-distribution
		(1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	}

	adapt <- function(tune,keep,k,target=0.44){  # adaptive tuning
		a <- min(0.01,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	get.Sigma <- function(sigma2,a,rho){  
		# Get covariance matrix, determinant, and precision matrix
		S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		list(P=P,b=b,S=S)
	}

	dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# browser()		
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
			# +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		if(!log) out <- 0.5*out*b^(-0.5)  # density
			# *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		out
	}

	# Test dmvt2 function
	# test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# test2 <- log(0.5)+
		# log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# plot(test1,test2)
	# summary(test2-test1)
	# lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# test2 <- 0.5*
		# (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# plot(test1,test2)
	# summary(test2/test1)
	# gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# Accept/reject proposals for mu
		# browser()
		mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
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
		mh.star <- num.z1+num.z0+U[mu.star[x],]%*%gamma  # numerator of Metropolis-Hastings ratio
			 #-log(sum(exp(U%*%gamma))))  # integral over S.tilde
			# log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		mh.0 <-	denom.z1+denom.z0+U[mu[x],]%*%gamma   # denominator of Metropolis-Hastings ratio
			#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		exp(mh.star-mh.0)>runif(1)  # Accept or reject
	}


	#####################################################################
	###  Setup Variables 
	#####################################################################
# browser() 
	cat("\nSetting up variables....")

	Ts <- nrow(s)  # number of telemetry locations
	Ty <- length(y)  # number of wet/dry observations
	qX <- ncol(X)
	qW <- ncol(W)
	qU <- dim(U)[3]+1	# number of raster layers plus an intercept
	v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	y1 <- which(y==1)+Ts
	y0 <- which(y==0)+Ts
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	W.cross <- t(W)%*%W  # cross product of W
	idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# information in S.tilde; note that mu below references row idx in S.tilde
	h <- match(start$h,S.tilde[,2])  # Note: h corresponds to row idx of S.tilde 
		# for computational efficiency, and not idx of mu as in Appendix A
	U <- cbind(1,values(U)[idx])  # convert raster to design matrix

	lc <- as.numeric(priors$lc)  # Argos location quality class
	n.lc <- length(unique(lc))  # number of error classes
	lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	
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
	priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu
	priors$mu.sigma <- priors$mu.sigma/s.sd  # variance on lognormal prior for sigma.mu

	# Center and scale starting values
	sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	sigma <- start$sigma/s.sd  # observation model standard deviation


	#####################################################################
	### Priors
	#####################################################################
	
	cat("\nGetting priors....")
	
	# Observation model
	u.sigma <- priors$u.sigma/s.sd
	
	# Temporal haul-out process model 
	mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# Spatial haul-out process model 
	J <- priors$J  # maximum number of clusters per truncation approximation
	mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	Sigma.gamma.inv <- solve(Sigma.gamma)


	#####################################################################
	### Appendix A, Step 1: starting values 
	#####################################################################

	cat("\nGetting starting values....")
# browser()
	# Observation model
	a <- start$a
	rho <- start$rho
	Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# Temporal haul-out process model
	beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)
  	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	z <- start$z  # haul-out status; z=1: dry, z=0: wet
	linpred <- X%*%beta+W%*%alpha

	# Spatial haul-out process model
	gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	theta <- start$theta  # DP concentration parameter


  	#####################################################################
	### Create receptacles for output
	#####################################################################
  
  	cat("\nCreating receptacles for output....")
	beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	theta.save <- numeric(n.mcmc)  # DP concentration parameter
	m.save <- numeric(n.mcmc)  # number of clusters
	v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

	sigma.save <- matrix(0,n.mcmc,n.lc)
	a.save <- matrix(0,n.mcmc,n.lc)
	rho.save <- matrix(0,n.mcmc,n.lc)

    
	#####################################################################
	### Appendix A, Steps 3-8: MCMC loop 
	#####################################################################

	keep <- list(mu=0,sigma.mu=0,gamma=0)  # number of proposals accepted for Metropolis updates

	keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))
	keep.tmp <- keep  # for adaptive tuning
	m.save.tmp <- 0
	t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	t.mcmc.start <- Sys.time()  # timing MCMC iterations
	T.b <- 50  # frequency of adaptive tuning
	
	cat("\nEntering MCMC Loop....\n")
	for (k in 1:n.mcmc) {
    	if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	cat(k,"");flush.console()	
    		plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	} 

    	if(k%%50==0) {  # Adaptive tuning
			# browser()
			keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			tune$sigma.mu <- adapt(tune$sigma.mu,keep.tmp$sigma.mu,k)
			tune$gamma <- adapt(tune$gamma,keep.tmp$gamma,k)
			tune$mu <- adapt(tune$mu,keep.tmp$mu,k)
			tune$sigma <- sapply(1:n.lc,function(x) adapt(tune$sigma[i],keep.tmp$sigma[i],k))
			tune$a <- sapply(1:n.lc,function(x) adapt(tune$a[i],keep.tmp$a[i],k))
			tune$rho <- sapply(1:n.lc,function(x) adapt(tune$rho[i],keep.tmp$rho[i],k))
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 
	
		#--------------------------------------------------------------------------
		# Appendix A, Step 4: update observation model parameters (Sigma)
		#--------------------------------------------------------------------------
# browser()
		sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		for(i in 1:n.lc){ #Loop to iterate over error classes: Appendix A, step 2(f)

			idx <- lc.list[[i]]  # index of locations in error class i
			z1 <- idx[which(z[idx]==1)]
			z0 <- idx[which(z[idx]==0)]

			### Sample sigma: Appendix A, step 2(b)

			if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				mh.star.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				mh.0.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					sigma[i] <- sigma.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$sigma[i] <- keep$sigma[i]+1
					keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				}
			}

			### Sample a: Appendix A, step 2(c)

			if(a.star[i]>0 & a.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				mh.star.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				mh.0.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				if(exp(mh.star.a-mh.0.a)>runif(1)){
					a[i] <- a.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$a[i] <- keep$a[i]+1
					keep.tmp$a[i] <- keep.tmp$a[i]+1
				}
			}

			### Sample rho: Appendix A, step 2(d)

			if(rho.star[i]>0 & rho.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				mh.star.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				mh.0.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
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
		# Appendix A, Step 4: update temporal haul-out process model parameters 
		#--------------------------------------------------------------------------
# browser()
	 	t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		###
		### Appendix A, Step 4(a): update v(t_y)
		###
		
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		###
		### Appendix A, Step 4(b): update v(t_s)
		###
			
		z1 <- which(z==1)
		z0 <- which(z==0)
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	
		###
		### Appendix A, Step 4(c): update alpha
		###
		
		b <- crossprod(W,(v-X%*%beta))
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		###	
		### Appendix A, Step 4(d): update beta
		###
		
		b <- crossprod(X,(v-W%*%alpha))
	  	beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of auxilliary variable update
	  	
	  	###
	  	### Appendix A, Step 4(e): calculate Prob(z(t_s)==1)
	  	###
	  	
	  	linpred <- X%*%beta+W%*%alpha  # update linear predictor 
		p <- pnorm(linpred[1:Ts,])

	    ###
	    ### Appendix A, Step 4(f): update z
	    ###
		
		p1 <- p*sapply(1:Ts,function(x) 
			dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		p2 <- (1-p)*sapply(1:Ts,function(x) 
			dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		p <- exp(log(p1)-log(p1+p2))	
		z <- rbinom(Ts,1,p)

		###
		### Appendix A, Step 4(g): update sigma2.alpha
		###
		
		r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		q.tmp <- qW/2+priors$q.sigma.alpha
		sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		diag(Sigma.alpha) <- sigma2.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		

	    #--------------------------------------------------------------------------
	  	# Appendix A, Step 5: update sigma.mu
	    #--------------------------------------------------------------------------
# browser()

		# Lognormal prior
	    sigma.mu.star <-  exp(rnorm(1,log(sigma.mu),tune$sigma.mu))
		Q.star <- sapply(1:n.lc,function(x) 
			get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
		idx <- which(z==0)
	    mh.star.sigma.mu <-	sum(sapply(idx,function(x) 
	    	dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
	    	dnorm(log(sigma.mu.star),log(priors$mu.sigma),priors$tau,log=TRUE)
	    mh.0.sigma.mu <- sum(sapply(idx,function(x)
	    	dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
	    	dnorm(log(sigma.mu),log(priors$mu.sigma),priors$tau,log=TRUE)
	    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	sigma.mu <- sigma.mu.star
			Q <- Q.star
			keep$sigma.mu <- keep$sigma.mu+1
			keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1

    	} 

		# Uniform prior	    
	    # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)				   
	    # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,n.lc,Mix)
			# Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,a,rho)
			# idx <- which(z==0)
		    # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		    # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Q <- Q.star
				# keep$sigma.mu <- keep$sigma.mu+1
	    	# } 
	    # }
		
				
		#--------------------------------------------------------------------------
		# Update spatial haul-out process model parameters: Appendix A, Step 6
		#--------------------------------------------------------------------------
# browser()
		###
		### Appendix A, Step 6(a): tabulate cluster membership 
		###

		n <- table(h)  
		m <- length(n)  # number of clusters
		mu <- as.numeric(names(n))  # idx of occupied clusters

		###
		### Appendix A, Step 6(d): update 'unoccupied' mu
		###

		p <- exp(U[-mu,]%*%gamma)
		mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		###
		### Appendix A, Step 6(c): update 'occupied' mu		   
		###
# browser()			
		mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		mh <- sapply(dup.idx,function(x)  # accepted proposals
			get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		keep$mu <- keep$mu+sum(mh)
		keep.tmp$mu <- keep.tmp$mu+sum(mh)
		keep.idx <- dup.idx[mh]
		mu[keep.idx] <- mu.star[keep.idx]  # update mu


		###
	    ### Appendix A, Step 6(f): update gamma
	    ###

		# Integral over S.tilde, occupied mu only
		gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))+
			sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
		mh.0.gamma <- sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))+
			sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    gamma <- gamma.star
	        # gamma.int <- gamma.int.star
	        keep$gamma <- keep$gamma+1
	        keep.tmp$gamma <- keep.tmp$gamma+1
    	} 

		###
		### Appendix A, Step 6(b): update the stick-breaking process 
		###
		
			# Appendix A, Step 6(b(i)): create index set I
			I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# Appendix A, Step 6(b(ii)): update eta
			n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# Appendix A, Step 6(b(ii)): update pi
		    pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # Appendix A, Step 6(b(iii)): update theta
			theta <- rgamma(1,priors$r.theta+J-1,priors$q.theta-sum(log(1-eta[-J])))  

			
		###
	    ### Appendix A, Step 6(e): update h(t_s)
	    ###

		mu <- c(mu[I],mu.tmp)
		h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))


		###
		###  Appendix A, Step 7: save samples 		   
		###
		
		sigma.save[k,] <- sigma*s.sd
		a.save[k,] <- a
		rho.save[k,] <- rho

		mu.save[,k] <- S.tilde[h,2]
		theta.save[k] <- theta    
		sigma.mu.save[k] <- sigma.mu
		sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		gamma.save[k,] <- gamma
		v.save[,k] <- v
		z.save[,k] <- z
		m.save[k] <- m
		m.save.tmp <- m.save.tmp+m
	}
  	
  	tune$sigma.mu <- tune$sigma.mu*s.sd
	tune$mu <- tune$mu*s.sd
  	sigma.mu.save <- sigma.mu.save*s.sd
  	t.mcmc.end <- Sys.time()

	#####################################################################
	### Write output
	#####################################################################
	  
	keep$sigma.mu <- keep$sigma.mu/n.mcmc
	keep$mu <- keep$mu/sum(m.save)
	keep$gamma <- keep$gamma/n.mcmc
	keep$sigma <- keep$sigma/n.mcmc
	keep$a <- keep$a/n.mcmc
	keep$rho <- keep$rho/n.mcmc

	cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 

	cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	cat(paste("\na acceptance rate:",round(keep$a,2))) 
	cat(paste("\nrho acceptance rate:",round(keep$rho,2))) 

	cat(paste("\n\nEnd time:",Sys.time()))
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	cat(paste("\nTime per MCMC iteration:",
		round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	list(beta=beta.save,gamma=gamma.save, alpha=alpha.save,
		mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		sigma.mu=sigma.mu.save, sigma.alpha=sigma.alpha.save,
		keep=keep,tune=tune,n.mcmc=n.mcmc,
		sigma=sigma.save,a=a.save,rho=rho.save)
}
