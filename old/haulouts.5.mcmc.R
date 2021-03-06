haulouts.5.mcmc <- function(s,lc,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	###
 	### Brian M. Brost (28 DEC 2015)
 	### Cleaned up version of haulouts.4.mcmc.R with estimation of telemetry error parameters
 	###
 	
 	###
 	### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	### y=ancillary data source containing binary wet/dry status;
 	### X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s and wet/dry status of ancillary data y;  
 	### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	### the location of haul-out sites; S.tilde=matrix summarizing support of haul-out sites 
 	###
 	
 	### See Appendix A of haul-outs manuscript for write-up of this model

	t.start <- Sys.time()
	cat(paste("Start time:",t.start,"\n"))

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

	adapt <- function(tune,keep,k,target=0.44){  # adaptive tuning
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
		mh.star <- num.z1+num.z0+U[mu.star[x],]%*%gamma  # numerator of Metropolis ratio
			 #-log(sum(exp(U%*%gamma))))  # integral over S.tilde
			# log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		mh.0 <-	denom.z1+denom.z0+U[mu[x],]%*%gamma   # denominator of Metropolis ratio
			#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		exp(mh.star-mh.0)>runif(1)  # Accept or reject
	}

	#####################################################################
	###  Get starting values from previous model
	#####################################################################

	k <- nrow(start$sigma)
	if(!is.null(k)){ 
		start$sigma <- start$sigma[k,]
		start$a <- start$a[k,]
		start$rho <- start$rho[k,]
		start$theta <- start$theta[k]+0.1				
		start$sigma.mu <- start$sigma.mu[k]
		start$gamma <- start$gamma[k,] 
		start$h <- start$mu[,k]
		start$h <- S.tilde[match(start$h,S.tilde[,2]),1]
		start$z <- start$z[,k]
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
	qU <- ncol(U)
	v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	y1 <- which(y==1)+Ts
	y0 <- which(y==0)+Ts
	y1.sum <- length(y1)
	y0.sum <- length(y0)
	W.cross <- t(W)%*%W  # cross product of W
	lc <- as.numeric(lc)  # Argos location quality class
	n.lc <- length(unique(lc))  # number of error classes
	lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	

	#####################################################################
	### Priors
	#####################################################################
	
	cat("\nGetting priors....")
	
	# Observation model
	u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	r.sigma.alpha <- priors$r.sigma.alpha  # IG prior on sigma.alpha
	q.sigma.alpha <- priors$q.sigma.alpha  # IG prior on sigma.alpha
	r.theta <- priors$r.theta  # IG prior on theta
	q.theta <- priors$q.theta  # IG prior on theta
	sigma.gamma <- priors$sigma.gamma  # variance of normal prior on gamma
	
	# Temporal haul-out process model 
	mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	Sigma.beta <- diag(qX)*priors$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)
	A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# Spatial haul-out process model 
	J <- priors$J  # maximum number of clusters per truncation approximation
	mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# Sigma.gamma <- diag(qU)*sigma.gamma^2
  	# Sigma.gamma.inv <- solve(Sigma.gamma)


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
	### Starting values: Appendix A, Steps 1 and 2
	#####################################################################

	cat("\nGetting starting values....")
# browser()
	# Observation model
	a <- start$a
	rho <- start$rho
	Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# Temporal haul-out process model
	alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	Sigma.alpha.inv <- solve(Sigma.alpha)
  	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	z <- start$z  # haul-out status; z=1: dry, z=0: wet

	# Spatial haul-out process model
	gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	theta <- start$theta  # DP concentration parameter


  	#####################################################################
	### Create receptacles for output
	#####################################################################
  
  	cat("\nCreating receptacles for output....")
	sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudianl telemetry measurement error
	a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
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
	
    
	#####################################################################
	### MCMC loop: Appendix A, Steps 3-9
	#####################################################################
	
	# Track overall MH accpetance rate
	keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# Track MH accpetance rate for adaptive tuning
	keep.tmp <- keep  
	T.b <- 50  # frequency of adaptive tuning
	m.save.tmp <- 0  # number of clusters

	t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	t.mcmc.start <- Sys.time()  # timing MCMC iterations
	
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
			tune$sigma <- sapply(1:n.lc,function(x) adapt(tune$sigma[x],keep.tmp$sigma[x],k))
			tune$a <- sapply(1:n.lc,function(x) adapt(tune$a[x],keep.tmp$a[x],k))
			tune$rho <- sapply(1:n.lc,function(x) adapt(tune$rho[x],keep.tmp$rho[x],k))
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
			m.save.tmp <- 0
	   	} 
	
	    
	    #--------------------------------------------------------------------------
	  	# Update sigma.mu: Appendix A, Step 4 
	    #--------------------------------------------------------------------------
# browser()
		# Lognormal prior
	    sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		Q.star <- sapply(1:n.lc,function(x) 
			get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
		idx <- which(z==0)
	    mh.star.sigma.mu <-	sum(sapply(idx,function(x) 
	    	dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
	    	dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
	    mh.0.sigma.mu <- sum(sapply(idx,function(x)
	    	dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
	    	dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
	    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	sigma.mu <- sigma.mu.star
			Q <- Q.star
			keep$sigma.mu <- keep$sigma.mu+1
			keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
    	} 

	
		#--------------------------------------------------------------------------
		# Update observation model parameters: Appendix A, Step 5
		#--------------------------------------------------------------------------
# browser()
		sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		for(i in 1:n.lc){  # loop to iterate over error classes: Appendix A, Step 5(b)

			idx <- lc.list[[i]]  # index of locations in error class i
			z1 <- idx[which(z[idx]==1)]
			z0 <- idx[which(z[idx]==0)]

			###
			### Sample sigma: Appendix A, Step 5(a.i)
			###

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

			###
			### Sample a: Appendix A, Step 5(a.ii)
			###
			
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

			###
			### Sample rho: Appendix A, Step 5(a.iii)
			###
			
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
		# Update temporal haul-out process model parameters: Appendix A, Step 6
		#--------------------------------------------------------------------------
# browser()
	 	t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		###	
		### Update beta: Appendix A, Step 6(a)
		###
		
		b <- crossprod(X,(v-W%*%alpha))
	  	beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		###
		### Update alpha: Appendix A, Step 6(b)
		###
		
		b <- crossprod(W,(v-X%*%beta))
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		###
		### Update v(t_y): Appendix A, Step 6(c)
		###

		linpred <- X%*%beta+W%*%alpha  # update linear predictor 
	  	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		###
		### Update v(t_s): Appendix A, Step 6(d)
		###
			
		z1 <- which(z==1)
		z0 <- which(z==0)
		v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of aux. variable update
	  	
	  	###
	  	### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	###
	  	
		p <- pnorm(linpred[1:Ts,])

	    ###
	    ### Update z: Appendix A, Step 6(f)
	    ###
		
		p1 <- p*sapply(1:Ts,function(x) 
			dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		p2 <- (1-p)*sapply(1:Ts,function(x) 
			dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		p <- exp(log(p1)-log(p1+p2))	
		z <- rbinom(Ts,1,p)

		###
		### Update sigma2.alpha: Appendix A, Step 6(g)
		###
		
		r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.sigma.alpha)
		q.tmp <- qW/2+q.sigma.alpha
		sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		diag(Sigma.alpha) <- sigma2.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		
					
		#--------------------------------------------------------------------------
		# Update spatial haul-out process model parameters: Appendix A, Step 7
		#--------------------------------------------------------------------------
# browser()
		###
		### Tabulate cluster membership: Appendix A, Step 7(a) 
		###

		n <- table(h)  
		m <- length(n)  # number of clusters
		mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		###
		### Update the stick-breaking process: Appendix A, Step 7(b)
		###
		
			# Create index set I: Appendix A, Step 6(b.i)
			I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# Update eta: Appendix A, Step 6(b.ii)
			n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# Update pi: Appendix A, Step 6(b.ii)
		    pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # Update theta: Appendix A, Step 6(b.iii)
			theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  

		###
		### Update 'unoccupied' mu: Appendix A, Step 7(c)
		###

		p <- exp(U[-mu,]%*%gamma)
		mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		###
		### Update 'occupied' mu: Appendix A, Step 7(d)		   
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
	    ### Update gamma: Appendix A, Step 7(e)
	    ###
	    
		# Integral over S.tilde, occupied mu only
		gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,sigma.gamma,log=TRUE))+
			sum(n*c(U[mu,]%*%gamma.star)-n*log(sum(exp(U%*%gamma.star))))
			# sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
		mh.0.gamma <- sum(dnorm(gamma,mu.gamma,sigma.gamma,log=TRUE))+
			sum(n*c(U[mu,]%*%gamma)-n*log(sum(exp(U%*%gamma))))
			# sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    gamma <- gamma.star
	        # gamma.int <- gamma.int.star
	        keep$gamma <- keep$gamma+1
	        keep.tmp$gamma <- keep.tmp$gamma+1
    	} 
			
		###
	    ### Update h(t_s): Appendix A, Step 7(f)
	    ###

		mu <- c(mu[I],mu.tmp)
		h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		#--------------------------------------------------------------------------
		# Save samples: Appendix A, Step 8
		#--------------------------------------------------------------------------
						
		sigma.save[k,] <- sigma
		a.save[k,] <- a
		rho.save[k,] <- rho
		mu.save[,k] <- S.tilde[h,2]
		theta.save[k] <- theta    
		sigma.mu.save[k] <- sigma.mu
		sigma.alpha.save[k] <- sigma2.alpha
		alpha.save[k,] <- alpha
		beta.save[k,] <- beta
		gamma.save[k,] <- gamma
		v.save[,k] <- v
		z.save[,k] <- z
		m.save[k] <- m
		m.save.tmp <- m.save.tmp+m
	}  # end of MCMC loop

  	t.mcmc.end <- Sys.time()
  	
  	tune$sigma.mu <- tune$sigma.mu*s.sd
	tune$mu <- tune$mu*s.sd
  	sigma.save <- sigma.save*s.sd
  	sigma.mu.save <- sigma.mu.save*s.sd
	sigma.alpha.save <- sqrt(sigma.alpha.save)
	
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
	cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	cat("\na acceptance rate:",round(keep$a,2))
	cat("\nrho acceptance rate:",round(keep$rho,2)) 

	cat(paste("\n\nEnd time:",Sys.time()))
	cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	cat(paste("\nTime per MCMC iteration:",
		round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	list(beta=beta.save,gamma=gamma.save,alpha=alpha.save,
		mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		sigma.mu=sigma.mu.save,sigma.alpha=sigma.alpha.save,
		sigma=sigma.save,a=a.save,rho=rho.save,
		keep=keep,tune=tune,n.mcmc=n.mcmc)
}
