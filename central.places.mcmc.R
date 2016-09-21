central.places.mcmc <- function(s,lc,y,X,W,J,S.tilde,priors,tune,start,
	adapt=FALSE,n.mcmc){
 
 	###
 	### Brian M. Brost (28 DEC 2015)
 	###
 	
 	###
 	### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	### y=ancillary data source containing binary wet/dry status;
 	### X=design matrix containing covariates influencing wet/dry
 	### status of telemetry locations s and wet/dry status of ancillary data y;  
 	### W=basis expansion for y; S.tilde=matrix summarizing support of haul-out sites;  
 	### priors=prior specifications; tune=tuning parameters;start=starting values;
 	### adapt=enable adaptive tuning; n.mcmc=number of MCMC iterations
 	###
 	
 	### See Appendix A of haul-outs manuscript for model statement and full-conditional
 	### distributions associated with this model, and Appendix B for corresponding
 	### pseudocode.

### This algorithm has been modified such that s and y are temporally aligned, and thus
### z=y. Modifications in the code are commented out and indicated by a "*", i.e., "#*"	

	t.start <- Sys.time()
	cat(paste("\n\nStart time:",t.start,"\n"))

	#####################################################################
	### Libraries and Subroutines
	#####################################################################
# browser()
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
			# +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		if(!log) out <- 0.5*out*b^(-0.5)  # density
			# *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
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
# browser() 
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

	#####################################################################
	### Priors
	#####################################################################
	
	cat("\nGetting priors....")
	
	# Observation model
	u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	r.alpha <- priors$r.alpha  # IG prior on sigma.alpha
	q.alpha <- priors$q.alpha  # IG prior on sigma.alpha
	r.theta <- priors$r.theta  # IG prior on theta
	q.theta <- priors$q.theta  # IG prior on theta
	
	# Temporal haul-out process model 
	mu.alpha <- matrix(0,qW,1)  # mean of random effects
	mu.beta <- matrix(0,qX,1)  # mean of coefficients on X (fixed effects)
	Sigma.beta <- diag(qX)*priors$sigma.beta^2  # variance for beta
	Sigma.beta.inv <- solve(Sigma.beta)
	A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# Spatial haul-out process model
	f.hat <- priors$f.hat  # prior on mu.j
	

	#####################################################################
	### Standardize parameters
	#####################################################################

	cat("\nStandardizing variables....")
# browser()	
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
# f.hat <- f.hat/s.sd  # prior on mu.j

	# Center and scale starting values
	sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	sigma <- start$sigma/s.sd  # observation model standard deviation


	#####################################################################
	### Starting values: Appendix A, Steps 1 and 2
	#####################################################################

	cat("\nGetting starting values....")

	# Observation model
	a <- start$a
	rho <- start$rho
	Sigma <- sapply(1:n.lc,function(x)
		get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	Q <- sapply(1:n.lc,function(x) 
		get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# Temporal haul-out process model
	alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	Sigma.alpha <- diag(qW)*start$sigma.alpha^2  # variance of random effects
  	Sigma.alpha.inv <- solve(Sigma.alpha)
	W.cross <- t(W)%*%W  # cross product of W
	# W.cross <- t(W[-(1:Ts)])%*%W[-(1:Ts)]  # use when estimating z
  	A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	
	# When telemetry locations and wet/dry data are temporally aligned:
	z <- y  # haul-out status of telemetry locations
	z1 <- which(z==1)
	z0 <- which(z==0)
	z1.lc <- tapply(z1,lc[z1],I)  # by location quality class
	z0.lc <- tapply(z0,lc[z0],I)  # by location quality class

lc <- as.numeric(lc)  # numerical Argos location quality class

	# Otherwise, get starting values for z and estimate below
	# Note: the matrices X and W should be stacked such that 
	# X=rbind(X(t_s),X(t_y)) and W=rbind(W(t_s),W(t_y))
	# z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# y1 <- y1+Ts
	# y0 <- y0+Ts
	# v <- numeric(Ts+Ty)

	# Spatial haul-out process model
# browser()
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
	# z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable	
    
    
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
		# Update spatial haul-out process model parameters: Appendix B, Step 4
		#--------------------------------------------------------------------------

		###
	    ### Update h(t): Appendix B, Step 4(a)
	    ###
# browser()
		# Indexes records in mu.j (1:J)
		h <- sapply(1:Ts,function(x) sample(1:J,1,prob= 
			exp(log(pi)+dt2(s[x,],S.tilde[mu.j,3:4],z[x],Sigma,Q,lc[x]))))

		###
		### Tabulate cluster membership: Appendix B, Step 4(b) 
		###

		n <- table(h)
		m <- length(n)  # number of clusters			
		j <- as.numeric(names(n))  # idx of occupied clusters
	
		###
		### Update the stick-breaking weights: Appendix B, Step 4(c)
		###

		n.tmp <- numeric(J)
		n.tmp[j] <- n
		eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights
		eta[-J] <- abs(eta[-J]-10^-10)  # fudge factor for preventing eta[1:(J-1)]=1 or 0

		###
		### Update the stick-breaking probabilities: Appendix B, Step 4(d)
		###
		
	    pi <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities
	
		###
		### Update the DP concentration parameter: Appendix B, Step 4(e)
		###

		theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
		
		# MH update for theta with uniform prior
		# theta.star <- rnorm(1,theta,0.1)
		# if(theta.star>0.05&theta.star<0.75){
			# mh.star.theta <- sum(dbeta(eta[-J],1,theta.star,log=TRUE))	
		    # mh.0.theta <- sum(dbeta(eta[-J],1,theta,log=TRUE))
		    # if(exp(mh.star.theta-mh.0.theta)>runif(1)){
	        	# theta <- theta.star
				# # keep$sigma.mu <- keep$sigma.mu+1
				# # keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	# } 
		# }
		
		###
		### Update 'occupied' mu: Appendix B, Step 4(f)		   
		###
# browser()
		# mu.j indexes rows in S.tilde
		mu.star[j] <- sapply(mu.j[j],
			function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*
			dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		idx <- j[which(!duplicated(mu.star[j]))]  # exclude duplicate proposals
		mh <- sapply(idx,function(x)  # accepted proposals
			get.mh.mu(x,mu.j,mu.star,S.tilde,h,s,lc,z,Sigma,Q,f.hat)) 
		keep$mu <- keep$mu+sum(mh)
		keep.tmp$mu <- keep.tmp$mu+sum(mh)
		keep.idx <- idx[mh]
		mu.j[keep.idx] <- mu.star[keep.idx]  # update mu

		###
		### Update 'unoccupied' mu: Appendix A, Step 4(g)
		###

		mu.j[-j] <- sample(S.tilde[-mu.j[j],1],J-m,replace=FALSE,prob=f.hat[-mu.j[j]])  
			# idx of unoccupied mu

		###
		### Use h to map mu.j to times t: Appendix B, Step 4(h)
		###

		mu.t <- S.tilde[mu.j[h],3:4]  # row idx of S.tilde


		#--------------------------------------------------------------------------
	  	# Update sigma.mu: Appendix B, Step 5 
	    #--------------------------------------------------------------------------
# browser()
		# Lognormal prior
# if(k>2000){	   
	    sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		if(sigma.mu.star>0){
			Q.star <- sapply(1:n.lc,function(x) 
				get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
			# z0 <- which(z==0)  # only necessary when estimating z
		    mh.star.sigma.mu <-	sum(sapply(z0,function(x) 
		    	dt2(s[x,],mu.t[x,],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    mh.0.sigma.mu <- sum(sapply(z0,function(x)
		    	dt2(s[x,],mu.t[x,],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	sigma.mu <- sigma.mu.star
				Q <- Q.star
				keep$sigma.mu <- keep$sigma.mu+1
				keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	} 
		}
# }
	
		#--------------------------------------------------------------------------
		# Update observation model parameters: Appendix B, Step 6
		#--------------------------------------------------------------------------
# browser()
		sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		for(i in 1:n.lc){  # iterate over error classes: Appendix B, Step 6(e)

			# Subset data for error class: Appendix B, Step 6(a)
			idx <- lc.list[[i]]  # index of locations in error class i
			# z1.tmp <- idx[which(z[idx]==1)]  # only necessary when estimating z
			# z0.tmp <- idx[which(z[idx]==0)]  # only necessary when estimating z
			z1.tmp <- z1.lc[[i]]
			z0.tmp <- z0.lc[[i]]
			s.z0 <- s[z0.tmp,]
			s.z1 <- s[z1.tmp,]
			mu.z0 <- mu.t[z0.tmp,]
			mu.z1 <- mu.t[z1.tmp,]

			###
			### Sample sigma: Appendix B, Step 6(b)
			###

			if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				mh.star.sigma <- sum(dt2(s.z1,mu.z1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma.star,Q.star))
				mh.0.sigma <- sum(dt2(s.z1,mu.z1,z=1,Sigma,Q,i))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma,Q,i))
				if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					sigma[i] <- sigma.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$sigma[i] <- keep$sigma[i]+1
					keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				}
			}

			###
			### Sample a: Appendix B, Step 6(c)
			###
			
			if(a.star[i]>0 & a.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				mh.star.a <- sum(dt2(s.z1,mu.z1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma.star,Q.star))
				mh.0.a <- sum(dt2(s.z1,mu.z1,z=1,Sigma,Q,i))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma,Q,i))
				if(exp(mh.star.a-mh.0.a)>runif(1)){
					a[i] <- a.star[i]
					Sigma[[i]] <- Sigma.star
					Q[[i]] <- Q.star
					keep$a[i] <- keep$a[i]+1
					keep.tmp$a[i] <- keep.tmp$a[i]+1
				}
			}

			###
			### Sample rho: Appendix B, Step 6(d)
			###
			
			if(rho.star[i]>0 & rho.star[i]<1){
				Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				mh.star.rho <- sum(dt2(s.z1,mu.z1,z=1,Sigma.star,Q.star))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma.star,Q.star))
				mh.0.rho <- sum(dt2(s.z1,mu.z1,z=1,Sigma,Q,i))+
					sum(dt2(s.z0,mu.z0,z=0,Sigma,Q,i))
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
		# Update temporal haul-out process model parameters: Appendix B, Step 7
		#--------------------------------------------------------------------------
# browser()
	 	t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		###	
		### Update beta: Appendix B, Step 7(a)
		###
		
		b <- crossprod(X,(v-W%*%alpha))
		# b <- crossprod(X[-(1:Ts),],(v[-(1:Ts)]-W[-(1:Ts),]%*%alpha))  # use when estimating z
		beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		###
		### Update alpha: Appendix B, Step 7(b)
		###
		
		b <- crossprod(W,(v-X%*%beta))
		# b <- crossprod(W[-(1:Ts),],(v[-(1:Ts)]-X[-(1:Ts),]%*%beta))  # use when estimating z
		alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		###
		### Update v(t): Appendix B, Step 7(c)
		###

		linpred <- X%*%beta+W%*%alpha  # update linear predictor 
	 	v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# End time of aux. variable update
	  	t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  

		###
		### Update sigma2.alpha: Appendix B, Step 7(d)
		###

		r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.alpha)
		q.tmp <- qW/2+q.alpha
		sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		diag(Sigma.alpha) <- sigma2.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)

		#--------------------------------------------------------------------------
		# Prediction of z(t_s): Appendix C
		#--------------------------------------------------------------------------
				
		# z1 <- which(z==1)  # only necessary when estimating z
		# z0 <- which(z==0)  # only necessary when estimating z
		# v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	# v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	###
	  	### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	###
	  	
		# p <- pnorm(v[1:Ts])

	    ###
	    ### Update z: Appendix A, Step 6(f)
	    ###
		
		# p1 <- p*sapply(1:Ts,function(x) 
			# dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		# p2 <- (1-p)*sapply(1:Ts,function(x) 
			# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		# p <- exp(log(p1)-log(p1+p2))	
		# z <- rbinom(Ts,1,p)

				
		#--------------------------------------------------------------------------
		# Update spatial haul-out process model parameters: Appendix A, Step 7
		#--------------------------------------------------------------------------

		# ###
		# ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# ###

		# n <- table(h)  
		# m <- length(n)  
		# mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		# ###
		# ### Update the stick-breaking process: Appendix A, Step 7(b)
		# ###
		
			# # Create index set I: Appendix A, Step 6(b.i)
			# I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # Update eta: Appendix A, Step 6(b.ii)
			# # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order

			# n.tmp <- c(n,rep(0,J-m-1))  # membership in decreasing order
			# eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # Update pi: Appendix A, Step 6(b.ii)
		    # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # Update theta: Appendix A, Step 6(b.iii)
			# theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
			

		# ###
		# ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# ###

		# mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE)  # idx of unoccupied mu
		
		# ###
		# ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# ###
# # browser()			
		# mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q)) 
		# keep$mu <- keep$mu+sum(mh)
		# keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# ###
	    # ### Update h(t_s): Appendix A, Step 7(f)
	    # ###

		# # mu <- c(mu[I],mu.tmp)
		# mu <- c(mu,mu.tmp)
		# h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

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
		# z.save[,k] <- z  # only necessary when estimating z
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

	# Update starting values
	start$sigma <- sigma*s.sd
	start$a <- a
	start$rho <- rho
	start$theta <- theta
	start$sigma.mu <- sigma.mu*s.sd
	start$mu <- S.tilde[mu.j,2]
	start$pi <- pi
	start$alpha <- alpha
	start$beta <- beta
	start$sigma.alpha <- sqrt(sigma2.alpha)
	# start$z <- start$z[,k]  # only necessary when estimating z
  
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
		v=v.save,beta=beta.save,alpha=alpha.save,  # z=z.save,
		sigma.alpha=sigma.alpha.save,sigma=sigma.save,a=a.save,rho=rho.save,
		keep=keep,tune=tune,priors,J=J,start=start,priors=priors,n.mcmc=n.mcmc)
}
















# haulouts.mcmc <- function(s,lc,y,X,W,J,S.tilde,priors,tune,start,
	# adapt=FALSE,n.mcmc){
 
 	# ###
 	# ### Brian M. Brost (28 DEC 2015)
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	# ### y=ancillary data source containing binary wet/dry status;
 	# ### X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# ### the location of haul-out sites; S.tilde=matrix summarizing support of haul-out sites 
 	# ###
 	
 	# ### See Appendix A of haul-outs manuscript for write-up of this model

	# ### This algorithm has been modified such that s and y are temporally aligned, and thus
	# ### z=y. Modifications in the code are commented out and indicated by a "*", i.e., "#*"	

	# t.start <- Sys.time()
	# cat(paste("\n\nStart time:",t.start,"\n"))

	# #####################################################################
	# ### Libraries and Subroutines
	# #####################################################################
# # browser()
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# }

	# get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# # a <- min(0.01,1/sqrt(k))
		# a <- min(0.025,1/sqrt(k))
		# exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# }

	# get.Sigma <- function(sigma2,a,rho){  # Get covariance matrix
		# S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		# b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		# P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		# list(P=P,b=b,S=S)
	# }

	# dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # Density of mixture t-distribution
		# x <- matrix(x,,2)	
		# y <- matrix(y,,2)
		# if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# if(!is.null(lc)){
			# S <- S[[lc]]
			# Q <- Q[[lc]]
		# }
		# P <- ifelse(z==1,S,Q)[[1]]  # precision matrix
		# b <- ifelse(z==1,S$b,Q$b)  # determinant
		# d <- x-y
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# out
	# }
  
	# get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q){
		# # Accept/reject proposals for mu
		# mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# idx.0 <- which(h==x&z==0)  # obs. associated with mu and z=0
		# idx.1 <- which(h==x&z==1)  # obs. associated with mu and z=1
		# num.z1 <- denom.z1 <- num.z0 <- denom.z0 <- 0
		# if(length(idx.1)>0){
			# num.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy.star,z=1,Sigma,Q,lc[x],log=TRUE)))
			# denom.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy,z=1,Sigma,Q,lc[x],log=TRUE)))
		# }
		# if(length(idx.0)>0){
			# num.z0 <- sum(sapply(idx.0,function(x) 
				# dt2(s[x,],mu.xy.star,z=0,Sigma,Q,lc[x],log=TRUE)))
			# denom.z0 <- sum(sapply(idx.0,function(x)
				# dt2(s[x,],mu.xy,z=0,Sigma,Q,lc[x],log=TRUE)))
		# }
		# mh.star <- num.z1+num.z0  # numerator of Metropolis ratio
		# mh.0 <-	denom.z1+denom.z0  # denominator of Metropolis ratio
		# exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# }

	# #####################################################################
	# ###  Get starting values from previous model
	# #####################################################################
# # browser()
	# k <- nrow(start$sigma)
	# if(!is.null(k)){ 
		# start$sigma <- start$sigma[k,]
		# start$a <- start$a[k,]
		# start$rho <- start$rho[k,]
		# start$theta <- start$theta[k]+0.1				
		# start$sigma.mu <- start$sigma.mu[k]
		# start$h <- start$mu[,k]
		# start$h <- S.tilde[match(start$h,S.tilde[,2]),1]
		# start$alpha <- start$alpha[k,]
		# start$beta <- start$beta[k,]
		# start$sigma.alpha <- start$sigma.alpha[k]
		# # start$z <- start$z[,k]  # only necessary when estimating z
	# }


	# #####################################################################
	# ###  Setup Variables 
	# #####################################################################
# # browser() 
	# cat("\nSetting up variables....")
	# Ts <- nrow(s)  # number of telemetry locations
	# Ty <- length(y)  # number of wet/dry observations
	# qX <- ncol(X)
	# qW <- ncol(W)
	# v <- numeric(Ty)  # auxilliary variable for haul-out process
	# y1 <- which(y==1)
	# y0 <- which(y==0)
	# y1.sum <- length(y1)
	# y0.sum <- length(y0)
	# lc <- as.numeric(lc)  # Argos location quality class
	# n.lc <- length(unique(lc))  # number of error classes
	# lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)


	# #####################################################################
	# ### Priors
	# #####################################################################
	
	# cat("\nGetting priors....")
	
	# # Observation model
	# u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	# mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	# sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	# r.sigma.alpha <- priors$r.sigma.alpha  # IG prior on sigma.alpha
	# q.sigma.alpha <- priors$q.sigma.alpha  # IG prior on sigma.alpha
	# r.theta <- priors$r.theta  # IG prior on theta
	# q.theta <- priors$q.theta  # IG prior on theta
	
	# # Temporal haul-out process model 
	# mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)


	# #####################################################################
	# ### Standardize parameters
	# #####################################################################

	# cat("\nStandardizing variables....")
	
	# # Center and scale s and S.tilde
	# s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.sd
	
	# # Center and scale tuning parameters
	# tune$sigma.mu <- tune$sigma.mu/s.sd
	# tune$mu <- tune$mu/s.sd
	# tune$sigma <- tune$sigma/s.sd

	# # Center and scale priors
	# mu.sigma <- mu.sigma/s.sd  
	# u.sigma <- u.sigma/s.sd

	# # Center and scale starting values
	# sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# sigma <- start$sigma/s.sd  # observation model standard deviation


	# #####################################################################
	# ### Starting values: Appendix A, Steps 1 and 2
	# #####################################################################

	# cat("\nGetting starting values....")

	# # Observation model
	# a <- start$a
	# rho <- start$rho
	# Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	# Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# # Temporal haul-out process model
	# alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	# Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# Sigma.alpha.inv <- solve(Sigma.alpha)
	# W.cross <- t(W)%*%W  # cross product of W
	# # W.cross <- t(W[-(1:Ts)])%*%W[-(1:Ts)]  # use when estimating z
  	# A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	
	# # When telemetry locations and wet/dry data are temporally aligned:
	# z <- y  # haul-out status of telemetry locations
	# z1 <- which(z==1)
	# z0 <- which(z==0)
	# z1.lc <- tapply(z1,lc[z1],I)  # by location quality class
	# z0.lc <- tapply(z0,lc[z0],I)  # by location quality class

	# # Otherwise, get starting values for z and estimate below
	# # Note: the matrices X and W should be stacked such that 
	# # X=rbind(X(t_s),X(t_y)) and W=rbind(W(t_s),W(t_y))
	# # z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # y1 <- y1+Ts
	# # y0 <- y0+Ts
	# # v <- numeric(Ts+Ty)

	# # Spatial haul-out process model
# # browser()
	# theta <- start$theta  # DP concentration parameter
	# h <- start$h
	# n <- table(h)
	# I <- order(n,decreasing=TRUE)
	# m <- length(n)
	# mu <- as.numeric(names(n[I]))
	# h.tmp <- match(h,mu)
	# mu <- c(mu,sample(S.tilde[,1],J-m,replace=FALSE))
	# n <- c(n,rep(0,J-m))
	# print(head(cbind(mu,n)))
	# # eta <- start$eta
	# # mu <- start$mu  # sample(S.tilde[,1],J,replace=FALSE)


  	# #####################################################################
	# ### Create receptacles for output
	# #####################################################################
  
  	# cat("\nCreating receptacles for output....")
	# sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudianl telemetry measurement error
	# a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	# rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
	# beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# m.save <- numeric(n.mcmc)  # number of clusters
	# v.save <- matrix(0,Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects
	# # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable	
    
    
	# #####################################################################
	# ### MCMC loop: Appendix A, Steps 3-9
	# #####################################################################
	
	# # Track overall MH accpetance rate
	# keep <- list(mu=0,sigma.mu=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# # Track MH accpetance rate for adaptive tuning
	# keep.tmp <- keep  
	# T.b <- 50  # frequency of adaptive tuning
	# m.save.tmp <- 0  # number of clusters

	# t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# t.mcmc.start <- Sys.time()  # timing MCMC iterations
	
	# cat("\nEntering MCMC Loop....\n")
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
	    	# print(head(cbind(mu,n)))
    	# } 
# # browser()

		# if(adapt==TRUE & k%%50==0) {  # Adaptive tuning
			# # browser()
			# keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# tune$sigma.mu <- get.tune(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# tune$mu <- get.tune(tune$mu,keep.tmp$mu,k)
			# tune$sigma <- sapply(1:n.lc,function(x) 
				# get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			# tune$a <- sapply(1:n.lc,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			# tune$rho <- sapply(1:n.lc,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			# keep.tmp <- lapply(keep.tmp,function(x) x*0)
			# m.save.tmp <- 0
	   	# } 	
		
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------
# # browser()

		# # # Update pi: Appendix A, Step 6(b.ii)
	    # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities
		
		# # ###
	    # # ### Update h(t): Appendix A, Step 7(f)
	    # # ###

		# # # Index records in mu
		# # h.tmp <- sapply(1:Ts,function(x) sample(1:J,1,prob= 
			# # exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# # ###
		# # ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# # ###
		
		# # n <- sapply(1:J,function(x) sum(h.tmp==x))  # tabulate cluster membership
		# # mu.idx <- which(n>0)  # idx of occupied clusters
		# # m <- sum(n>0)  # number of clusters

		# # # Update eta: Appendix A, Step 6(b.ii)
		# # eta <- c(rbeta(J-1,1+n,theta+Ts-cumsum(n)),1)  # stick-breaking weights

	    # # # Update theta: Appendix A, Step 6(b.iii)
		# # theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
		
		# # ###
		# # ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# # ###
# # # browser()			
		# # mu.star <- numeric(J)
		# # mu.star[mu.idx] <- sapply(mu[mu.idx],
			# # function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  

		# # dup.idx <- which((mu!=mu.star)&!duplicated(mu.star))  # exclude duplicate proposals
		# # dup.idx <- dup.idx[dup.idx%in%mu.idx]
		# # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h.tmp,s,lc,z,Sigma,Q)) 
		# # keep$mu <- keep$mu+sum(mh)
		# # keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# # keep.idx <- dup.idx[mh]
		# # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # ###
		# # ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# # ###

		# # mu[-mu.idx] <- sample(S.tilde[-mu[mu.idx],1],J-m,replace=FALSE)  # idx of unoccupied mu

		# # ###
		# # ### Map h (idx of mu) to mu (row idx of S.tilde)
		# # ###
		
		# # h <- mu[h.tmp]

	    
	    # #--------------------------------------------------------------------------
	  	# # Update sigma.mu: Appendix A, Step 4 
	    # #--------------------------------------------------------------------------
# # browser()
		# # Lognormal prior
	    # sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		# if(sigma.mu.star>0){
			# Q.star <- sapply(1:n.lc,function(x) 
				# get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
			# # z0 <- which(z==0)  # only necessary when estimating z
		    # mh.star.sigma.mu <-	sum(sapply(z0,function(x) 
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    # mh.0.sigma.mu <- sum(sapply(z0,function(x)
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Q <- Q.star
				# keep$sigma.mu <- keep$sigma.mu+1
				# keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	# } 
		# }

	
		# #--------------------------------------------------------------------------
		# # Update observation model parameters: Appendix A, Step 5
		# #--------------------------------------------------------------------------
# # browser()
		# sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		# a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		# rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		# for(i in 1:n.lc){  # loop to iterate over error classes: Appendix A, Step 5(b)

			# idx <- lc.list[[i]]  # index of locations in error class i
			# # z1.tmp <- idx[which(z[idx]==1)]  # only necessary when estimating z
			# # z0.tmp <- idx[which(z[idx]==0)]  # only necessary when estimating z
			# z1.tmp <- z1.lc[[i]]
			# z0.tmp <- z0.lc[[i]]

			# ###
			# ### Sample sigma: Appendix A, Step 5(a.i)
			# ###

			# if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				# Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# mh.star.sigma <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.sigma <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# sigma[i] <- sigma.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$sigma[i] <- keep$sigma[i]+1
					# keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				# }
			# }

			# ###
			# ### Sample a: Appendix A, Step 5(a.ii)
			# ###
			
			# if(a.star[i]>0 & a.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# mh.star.a <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.a <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.a-mh.0.a)>runif(1)){
					# a[i] <- a.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$a[i] <- keep$a[i]+1
					# keep.tmp$a[i] <- keep.tmp$a[i]+1
				# }
			# }

			# ###
			# ### Sample rho: Appendix A, Step 5(a.iii)
			# ###
			
			# if(rho.star[i]>0 & rho.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# mh.star.rho <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.rho <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# rho[i] <- rho.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$rho[i] <- keep$rho[i]+1
					# keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				# }
			# }
		# }

	
		# #--------------------------------------------------------------------------
		# # Update temporal haul-out process model parameters: Appendix A, Step 6
		# #--------------------------------------------------------------------------
# # browser()
	 	# t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# ###	
		# ### Update beta: Appendix A, Step 6(a)
		# ###
		
		# b <- crossprod(X,(v-W%*%alpha))
		# # b <- crossprod(X[-(1:Ts),],(v[-(1:Ts)]-W[-(1:Ts),]%*%alpha))  # use when estimating z
		# beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		# ###
		# ### Update alpha: Appendix A, Step 6(b)
		# ###
		
		# b <- crossprod(W,(v-X%*%beta))
		# # b <- crossprod(W[-(1:Ts),],(v[-(1:Ts)]-X[-(1:Ts),]%*%beta))  # use when estimating z
		# alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# ###
		# ### Update v(t): Appendix A, Step 6(c)
		# ###

		# linpred <- X%*%beta+W%*%alpha  # update linear predictor 
	 	# v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	# v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# # End time of aux. variable update
	  	# t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  

		# ###
		# ### Update sigma2.alpha: Appendix A, Step 6(g)
		# ###

		# r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.sigma.alpha)
		# q.tmp <- qW/2+q.sigma.alpha
		# sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# diag(Sigma.alpha) <- sigma2.alpha
		# Sigma.alpha.inv <- solve(Sigma.alpha)
		# A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)

		# ###
		# ### Prediction of z(t_s): Appendix C
		# ###
				
		# # z1 <- which(z==1)  # only necessary when estimating z
		# # z0 <- which(z==0)  # only necessary when estimating z
		# # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	# # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	# ###
	  	# ### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	# ###
	  	
		# # p <- pnorm(v[1:Ts])

	    # ###
	    # ### Update z: Appendix A, Step 6(f)
	    # ###
		
		# # p1 <- p*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		# # p2 <- (1-p)*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		# # p <- exp(log(p1)-log(p1+p2))	
		# # z <- rbinom(Ts,1,p)


		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------

		# ###
		# ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# ###

		# n <- sapply(1:J,function(x) sum(h.tmp==x))  # tabulate cluster membership
		# mu.idx <- which(n>0)  # idx of occupied clusters
		# m <- sum(n>0)  # number of clusters
		
		# ###
		# ### Update the stick-breaking process: Appendix A, Step 7(b)
		# ###

		# # Update eta: Appendix A, Step 6(b.ii)
		# eta <- c(rbeta(J-1,1+n,theta+Ts-cumsum(n)),1)  # stick-breaking weights

		# # Update pi: Appendix A, Step 6(b.ii)
	    # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities
		
		# # Update theta: Appendix A, Step 6(b.iii)
		# theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
	
		# ###
		# ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# ###
# # browser()			
		# mu.star <- numeric(J)
		# mu.star[mu.idx] <- sapply(mu[mu.idx],
			# function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# dup.idx <- which((mu!=mu.star)&!duplicated(mu.star))  # exclude duplicate proposals
		# dup.idx <- dup.idx[dup.idx%in%mu.idx]
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h.tmp,s,lc,z,Sigma,Q)) 
		# keep$mu <- keep$mu+sum(mh)
		# keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# ###
		# ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# ###

		# mu[-mu.idx] <- sample(S.tilde[-mu[mu.idx],1],J-m,replace=FALSE)  # idx of unoccupied mu

		# ###
	    # ### Update h(t): Appendix A, Step 7(f)
	    # ###

		# # Index records in mu
		# h.tmp <- sapply(1:Ts,function(x) sample(1:J,1,prob= 
			# exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# ###
		# ### Map h (idx of mu) to mu (row idx of S.tilde)
		# ###
		
		# h <- mu[h.tmp]

				
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------

		# # ###
		# # ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# # ###

		# # n <- table(h)  
		# # m <- length(n)  
		# # mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		# # ###
		# # ### Update the stick-breaking process: Appendix A, Step 7(b)
		# # ###
		
			# # # Create index set I: Appendix A, Step 6(b.i)
			# # I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # # Update eta: Appendix A, Step 6(b.ii)
			# # # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order

			# # n.tmp <- c(n,rep(0,J-m-1))  # membership in decreasing order
			# # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # Update pi: Appendix A, Step 6(b.ii)
		    # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # Update theta: Appendix A, Step 6(b.iii)
			# # theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
			

		# # ###
		# # ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# # ###

		# # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE)  # idx of unoccupied mu
		
		# # ###
		# # ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# # ###
# # # browser()			
		# # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q)) 
		# # keep$mu <- keep$mu+sum(mh)
		# # keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# # keep.idx <- dup.idx[mh]
		# # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # ###
	    # # ### Update h(t_s): Appendix A, Step 7(f)
	    # # ###

		# # # mu <- c(mu[I],mu.tmp)
		# # mu <- c(mu,mu.tmp)
		# # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# #--------------------------------------------------------------------------
		# # Save samples: Appendix A, Step 8
		# #--------------------------------------------------------------------------
						
		# sigma.save[k,] <- sigma
		# a.save[k,] <- a
		# rho.save[k,] <- rho
		# mu.save[,k] <- S.tilde[h,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# sigma.alpha.save[k] <- sigma2.alpha
		# alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# v.save[,k] <- v
		# m.save[k] <- m
		# m.save.tmp <- m.save.tmp+m
		# # z.save[,k] <- z  # only necessary when estimating z
	# }  # End of MCMC loop

  	# t.mcmc.end <- Sys.time()
  	
  	# tune$sigma.mu <- tune$sigma.mu*s.sd
  	# tune$sigma <- tune$sigma*s.sd
	# tune$mu <- tune$mu*s.sd
  	# sigma.save <- sigma.save*s.sd
  	# sigma.mu.save <- sigma.mu.save*s.sd
	# sigma.alpha.save <- sqrt(sigma.alpha.save)
	
  	# #####################################################################
	# ### Write output
	# #####################################################################
	  
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu <- keep$mu/sum(m.save)
	# keep$sigma <- keep$sigma/n.mcmc
	# keep$a <- keep$a/n.mcmc
	# keep$rho <- keep$rho/n.mcmc

	# cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	# cat("\na acceptance rate:",round(keep$a,2))
	# cat("\nrho acceptance rate:",round(keep$rho,2)) 

	# cat(paste("\n\nEnd time:",Sys.time()))
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# cat(paste("\nTime per MCMC iteration:",
		# round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	# list(mu=mu.save,theta=theta.save,m=m.save,sigma.mu=sigma.mu.save,
		# v=v.save,beta=beta.save,alpha=alpha.save,  # z=z.save,
		# sigma.alpha=sigma.alpha.save,sigma=sigma.save,a=a.save,rho=rho.save,
		# keep=keep,tune=tune,priors,J=J,n.mcmc=n.mcmc)
# }















# haulouts.mcmc <- function(s,lc,y,X,W,J,S.tilde,priors,tune,start,
	# adapt=FALSE,n.mcmc){
 
 	# ###
 	# ### Brian M. Brost (28 DEC 2015)
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	# ### y=ancillary data source containing binary wet/dry status;
 	# ### X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# ### the location of haul-out sites; S.tilde=matrix summarizing support of haul-out sites 
 	# ###
 	
 	# ### See Appendix A of haul-outs manuscript for write-up of this model

	# ### This algorithm has been modified such that s and y are temporally aligned, and thus
	# ### z=y. Modifications in the code are commented out and indicated by a "*", i.e., "#*"	

	# t.start <- Sys.time()
	# cat(paste("\n\nStart time:",t.start,"\n"))

	# #####################################################################
	# ### Libraries and Subroutines
	# #####################################################################
# # browser()
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# }

	# get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# # a <- min(0.01,1/sqrt(k))
		# a <- min(0.025,1/sqrt(k))
		# exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# }

	# get.Sigma <- function(sigma2,a,rho){  # Get covariance matrix
		# S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		# b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		# P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		# list(P=P,b=b,S=S)
	# }

	# dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # Density of mixture t-distribution
		# x <- matrix(x,,2)	
		# y <- matrix(y,,2)
		# if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# if(!is.null(lc)){
			# S <- S[[lc]]
			# Q <- Q[[lc]]
		# }
		# P <- ifelse(z==1,S,Q)[[1]]  # precision matrix
		# b <- ifelse(z==1,S$b,Q$b)  # determinant
		# d <- x-y
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# out
	# }
  
	# get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q){
		# # Accept/reject proposals for mu
		# mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# idx.0 <- which(h==x&z==0)  # obs. associated with mu and z=0
		# idx.1 <- which(h==x&z==1)  # obs. associated with mu and z=1
		# num.z1 <- denom.z1 <- num.z0 <- denom.z0 <- 0
		# if(length(idx.1)>0){
			# num.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy.star,z=1,Sigma,Q,lc[x],log=TRUE)))
			# denom.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy,z=1,Sigma,Q,lc[x],log=TRUE)))
		# }
		# if(length(idx.0)>0){
			# num.z0 <- sum(sapply(idx.0,function(x) 
				# dt2(s[x,],mu.xy.star,z=0,Sigma,Q,lc[x],log=TRUE)))
			# denom.z0 <- sum(sapply(idx.0,function(x)
				# dt2(s[x,],mu.xy,z=0,Sigma,Q,lc[x],log=TRUE)))
		# }
		# mh.star <- num.z1+num.z0  # numerator of Metropolis ratio
		# mh.0 <-	denom.z1+denom.z0  # denominator of Metropolis ratio
		# exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# }

	# #####################################################################
	# ###  Get starting values from previous model
	# #####################################################################
# browser()
	# k <- nrow(start$sigma)
	# if(!is.null(k)){ 
		# start$sigma <- start$sigma[k,]
		# start$a <- start$a[k,]
		# start$rho <- start$rho[k,]
		# start$theta <- start$theta[k]+0.1				
		# start$sigma.mu <- start$sigma.mu[k]
		# start$h <- start$mu[,k]
		# start$h <- S.tilde[match(start$h,S.tilde[,2]),1]
		# start$alpha <- start$alpha[k,]
		# start$beta <- start$beta[k,]
		# start$sigma.alpha <- start$sigma.alpha[k]
		# # start$z <- start$z[,k]  # only necessary when estimating z
	# }


	# #####################################################################
	# ###  Setup Variables 
	# #####################################################################
# # browser() 
	# cat("\nSetting up variables....")
	# Ts <- nrow(s)  # number of telemetry locations
	# Ty <- length(y)  # number of wet/dry observations
	# qX <- ncol(X)
	# qW <- ncol(W)
	# v <- numeric(Ty)  # auxilliary variable for haul-out process
	# y1 <- which(y==1)
	# y0 <- which(y==0)
	# y1.sum <- length(y1)
	# y0.sum <- length(y0)
	# lc <- as.numeric(lc)  # Argos location quality class
	# n.lc <- length(unique(lc))  # number of error classes
	# lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)


	# #####################################################################
	# ### Priors
	# #####################################################################
	
	# cat("\nGetting priors....")
	
	# # Observation model
	# u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	# mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	# sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	# r.sigma.alpha <- priors$r.sigma.alpha  # IG prior on sigma.alpha
	# q.sigma.alpha <- priors$q.sigma.alpha  # IG prior on sigma.alpha
	# r.theta <- priors$r.theta  # IG prior on theta
	# q.theta <- priors$q.theta  # IG prior on theta
	
	# # Temporal haul-out process model 
	# mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)


	# #####################################################################
	# ### Standardize parameters
	# #####################################################################

	# cat("\nStandardizing variables....")
	
	# # Center and scale s and S.tilde
	# s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.sd
	
	# # Center and scale tuning parameters
	# tune$sigma.mu <- tune$sigma.mu/s.sd
	# tune$mu <- tune$mu/s.sd
	# tune$sigma <- tune$sigma/s.sd

	# # Center and scale priors
	# mu.sigma <- mu.sigma/s.sd  
	# u.sigma <- u.sigma/s.sd

	# # Center and scale starting values
	# sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# sigma <- start$sigma/s.sd  # observation model standard deviation


	# #####################################################################
	# ### Starting values: Appendix A, Steps 1 and 2
	# #####################################################################

	# cat("\nGetting starting values....")

	# # Observation model
	# a <- start$a
	# rho <- start$rho
	# Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	# Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# # Temporal haul-out process model
	# alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	# Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# Sigma.alpha.inv <- solve(Sigma.alpha)
	# W.cross <- t(W)%*%W  # cross product of W
	# # W.cross <- t(W[-(1:Ts)])%*%W[-(1:Ts)]  # use when estimating z
  	# A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	
	# # When telemetry locations and wet/dry data are temporally aligned:
	# z <- y  # haul-out status of telemetry locations
	# z1 <- which(z==1)
	# z0 <- which(z==0)
	# z1.lc <- tapply(z1,lc[z1],I)  # by location quality class
	# z0.lc <- tapply(z0,lc[z0],I)  # by location quality class

	# # Otherwise, get starting values for z and estimate below
	# # Note: the matrices X and W should be stacked such that 
	# # X=rbind(X(t_s),X(t_y)) and W=rbind(W(t_s),W(t_y))
	# # z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # y1 <- y1+Ts
	# # y0 <- y0+Ts
	# # v <- numeric(Ts+Ty)

	# # Spatial haul-out process model
	# theta <- start$theta  # DP concentration parameter
	# h <- start$h

	# # eta <- start$eta
	# # mu <- start$mu  # sample(S.tilde[,1],J,replace=FALSE)


  	# #####################################################################
	# ### Create receptacles for output
	# #####################################################################
  
  	# cat("\nCreating receptacles for output....")
	# sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudianl telemetry measurement error
	# a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	# rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
	# beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# m.save <- numeric(n.mcmc)  # number of clusters
	# v.save <- matrix(0,Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects
	# # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable	
    
    
	# #####################################################################
	# ### MCMC loop: Appendix A, Steps 3-9
	# #####################################################################
	
	# # Track overall MH accpetance rate
	# keep <- list(mu=0,sigma.mu=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# # Track MH accpetance rate for adaptive tuning
	# keep.tmp <- keep  
	# T.b <- 50  # frequency of adaptive tuning
	# m.save.tmp <- 0  # number of clusters

	# t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# t.mcmc.start <- Sys.time()  # timing MCMC iterations
	
	# cat("\nEntering MCMC Loop....\n")
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# } 
# # browser()

		# if(adapt==TRUE & k%%50==0) {  # Adaptive tuning
			# # browser()
			# keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# tune$sigma.mu <- get.tune(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# tune$mu <- get.tune(tune$mu,keep.tmp$mu,k)
			# tune$sigma <- sapply(1:n.lc,function(x) 
				# get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			# tune$a <- sapply(1:n.lc,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			# tune$rho <- sapply(1:n.lc,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			# keep.tmp <- lapply(keep.tmp,function(x) x*0)
			# m.save.tmp <- 0
	   	# } 	
		
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------
# # browser()

		# # Update pi: Appendix A, Step 6(b.ii)
	    # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

# # pie[45:50] <- 0
		
		# ###
	    # ### Update h(t): Appendix A, Step 7(f)
	    # ###

		# h <- sapply(1:Ts,function(x) sample(1:J,1,prob= 
			# exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

# # table(h)

		# ###
		# ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# ###
		
		# n <- sapply(1:J,function(x) sum(h==x))  # tabulate cluster membership
		# mu.idx <- which(n>0)  # idx of occupied clusters
		# m <- sum(n>0)  # number of clusters

		# # Update eta: Appendix A, Step 6(b.ii)
		# eta <- c(rbeta(J-1,1+n,theta+Ts-cumsum(n)),1)  # stick-breaking weights

	    # # Update theta: Appendix A, Step 6(b.iii)
		# theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
		
		# ###
		# ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# ###
# # browser()			

		# mu.star <- numeric(J)
		# mu.star[mu.idx] <- sapply(mu[mu.idx],
			# function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  

		# dup.idx <- which((mu!=mu.star)&!duplicated(mu.star))  # exclude duplicate proposals
		# dup.idx <- dup.idx[dup.idx%in%mu.idx]
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q)) 
		# keep$mu <- keep$mu+sum(mh)
		# keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# ###
		# ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# ###

		# mu[-mu.idx] <- sample(S.tilde[-mu[mu.idx],1],J-m,replace=FALSE)  # idx of unoccupied mu

		
		# ###
		# ### Map h (idx of mu) to mu (row idx of S.tilde)
		# ###
		
		# h <- mu[h]

	    
	    # #--------------------------------------------------------------------------
	  	# # Update sigma.mu: Appendix A, Step 4 
	    # #--------------------------------------------------------------------------
# # browser()
		# # Lognormal prior
	    # sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		# if(sigma.mu.star>0){
			# Q.star <- sapply(1:n.lc,function(x) 
				# get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
			# # z0 <- which(z==0)  # only necessary when estimating z
		    # mh.star.sigma.mu <-	sum(sapply(z0,function(x) 
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    # mh.0.sigma.mu <- sum(sapply(z0,function(x)
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Q <- Q.star
				# keep$sigma.mu <- keep$sigma.mu+1
				# keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	# } 
		# }

	
		# #--------------------------------------------------------------------------
		# # Update observation model parameters: Appendix A, Step 5
		# #--------------------------------------------------------------------------
# # browser()
		# sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		# a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		# rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		# for(i in 1:n.lc){  # loop to iterate over error classes: Appendix A, Step 5(b)

			# idx <- lc.list[[i]]  # index of locations in error class i
			# # z1.tmp <- idx[which(z[idx]==1)]  # only necessary when estimating z
			# # z0.tmp <- idx[which(z[idx]==0)]  # only necessary when estimating z
			# z1.tmp <- z1.lc[[i]]
			# z0.tmp <- z0.lc[[i]]

			# ###
			# ### Sample sigma: Appendix A, Step 5(a.i)
			# ###

			# if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				# Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# mh.star.sigma <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.sigma <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# sigma[i] <- sigma.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$sigma[i] <- keep$sigma[i]+1
					# keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				# }
			# }

			# ###
			# ### Sample a: Appendix A, Step 5(a.ii)
			# ###
			
			# if(a.star[i]>0 & a.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# mh.star.a <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.a <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.a-mh.0.a)>runif(1)){
					# a[i] <- a.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$a[i] <- keep$a[i]+1
					# keep.tmp$a[i] <- keep.tmp$a[i]+1
				# }
			# }

			# ###
			# ### Sample rho: Appendix A, Step 5(a.iii)
			# ###
			
			# if(rho.star[i]>0 & rho.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# mh.star.rho <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma.star,Q.star))
				# mh.0.rho <- sum(dt2(s[z1.tmp,],S.tilde[h[z1.tmp],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0.tmp,],S.tilde[h[z0.tmp],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# rho[i] <- rho.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$rho[i] <- keep$rho[i]+1
					# keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				# }
			# }
		# }

	
		# #--------------------------------------------------------------------------
		# # Update temporal haul-out process model parameters: Appendix A, Step 6
		# #--------------------------------------------------------------------------
# # browser()
	 	# t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# ###	
		# ### Update beta: Appendix A, Step 6(a)
		# ###
		
		# b <- crossprod(X,(v-W%*%alpha))
		# # b <- crossprod(X[-(1:Ts),],(v[-(1:Ts)]-W[-(1:Ts),]%*%alpha))  # use when estimating z
		# beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		# ###
		# ### Update alpha: Appendix A, Step 6(b)
		# ###
		
		# b <- crossprod(W,(v-X%*%beta))
		# # b <- crossprod(W[-(1:Ts),],(v[-(1:Ts)]-X[-(1:Ts),]%*%beta))  # use when estimating z
		# alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# ###
		# ### Update v(t): Appendix A, Step 6(c)
		# ###

		# linpred <- X%*%beta+W%*%alpha  # update linear predictor 
	 	# v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	# v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# # End time of aux. variable update
	  	# t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  

		# ###
		# ### Update sigma2.alpha: Appendix A, Step 6(g)
		# ###

		# r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.sigma.alpha)
		# q.tmp <- qW/2+q.sigma.alpha
		# sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# diag(Sigma.alpha) <- sigma2.alpha
		# Sigma.alpha.inv <- solve(Sigma.alpha)
		# A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)

		# ###
		# ### Prediction of z(t_s): Appendix C
		# ###
				
		# # z1 <- which(z==1)  # only necessary when estimating z
		# # z0 <- which(z==0)  # only necessary when estimating z
		# # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	# # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	# ###
	  	# ### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	# ###
	  	
		# # p <- pnorm(v[1:Ts])

	    # ###
	    # ### Update z: Appendix A, Step 6(f)
	    # ###
		
		# # p1 <- p*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		# # p2 <- (1-p)*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		# # p <- exp(log(p1)-log(p1+p2))	
		# # z <- rbinom(Ts,1,p)

				
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------

		# # ###
		# # ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# # ###

		# # n <- table(h)  
		# # m <- length(n)  
		# # mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		# # ###
		# # ### Update the stick-breaking process: Appendix A, Step 7(b)
		# # ###
		
			# # # Create index set I: Appendix A, Step 6(b.i)
			# # I <- order(n,decreasing=TRUE)  # clusters ordered by membership

# # # browser()
# # # library(Hmisc)
# # # test <- rMultinom(rbind(c(.1,.2,.3,.4),c(.4,.3,.2,.1)),1)
# # # table(test)
		
			# # # Update eta: Appendix A, Step 6(b.ii)
			# # # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order

			# # n.tmp <- c(n,rep(0,J-m-1))  # membership in decreasing order
			# # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # Update pi: Appendix A, Step 6(b.ii)
		    # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # Update theta: Appendix A, Step 6(b.iii)
			# # theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  
			

		# # ###
		# # ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# # ###

		# # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE)  # idx of unoccupied mu
		
		# # ###
		# # ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# # ###
# # # browser()			
		# # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q)) 
		# # keep$mu <- keep$mu+sum(mh)
		# # keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# # keep.idx <- dup.idx[mh]
		# # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # ###
	    # # ### Update h(t_s): Appendix A, Step 7(f)
	    # # ###

		# # # mu <- c(mu[I],mu.tmp)
		# # mu <- c(mu,mu.tmp)
		# # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# #--------------------------------------------------------------------------
		# # Save samples: Appendix A, Step 8
		# #--------------------------------------------------------------------------
						
		# sigma.save[k,] <- sigma
		# a.save[k,] <- a
		# rho.save[k,] <- rho
		# mu.save[,k] <- S.tilde[h,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# sigma.alpha.save[k] <- sigma2.alpha
		# alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# v.save[,k] <- v
		# m.save[k] <- m
		# m.save.tmp <- m.save.tmp+m
		# # z.save[,k] <- z  # only necessary when estimating z
	# }  # End of MCMC loop

  	# t.mcmc.end <- Sys.time()
  	
  	# tune$sigma.mu <- tune$sigma.mu*s.sd
  	# tune$sigma <- tune$sigma*s.sd
	# tune$mu <- tune$mu*s.sd
  	# sigma.save <- sigma.save*s.sd
  	# sigma.mu.save <- sigma.mu.save*s.sd
	# sigma.alpha.save <- sqrt(sigma.alpha.save)
	
  	# #####################################################################
	# ### Write output
	# #####################################################################
	  
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu <- keep$mu/sum(m.save)
	# keep$sigma <- keep$sigma/n.mcmc
	# keep$a <- keep$a/n.mcmc
	# keep$rho <- keep$rho/n.mcmc

	# cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	# cat("\na acceptance rate:",round(keep$a,2))
	# cat("\nrho acceptance rate:",round(keep$rho,2)) 

	# cat(paste("\n\nEnd time:",Sys.time()))
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# cat(paste("\nTime per MCMC iteration:",
		# round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	# list(mu=mu.save,theta=theta.save,m=m.save,sigma.mu=sigma.mu.save,
		# v=v.save,beta=beta.save,alpha=alpha.save,  # z=z.save,
		# sigma.alpha=sigma.alpha.save,sigma=sigma.save,a=a.save,rho=rho.save,
		# keep=keep,tune=tune,priors,J=J,n.mcmc=n.mcmc)
# }






###
### MCMC algorithm for estimation of gamma, too
###

# haulouts.mcmc <- function(s,lc,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,
	# adapt=FALSE,n.mcmc,n.cores=NULL){
 
 	# ###
 	# ### Brian M. Brost (28 DEC 2015)
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	# ### y=ancillary data source containing binary wet/dry status;
 	# ### X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# ### the location of haul-out sites; S.tilde=matrix summarizing support of haul-out sites 
 	# ###
 	
 	# ### See Appendix A of haul-outs manuscript for write-up of this model

	# ### This algorithm has been modified such that s and y are temporally aligned, and thus
	# ### z=y. Modifications in the code are commented out and indicated by a "*", i.e., "#*"	

	# t.start <- Sys.time()
	# cat(paste("\n\nStart time:",t.start,"\n"))

	# #####################################################################
	# ### Libraries and Subroutines
	# #####################################################################
# # browser()
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# }

	# get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# # a <- min(0.01,1/sqrt(k))
		# a <- min(0.025,1/sqrt(k))
		# exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# }

	# get.Sigma <- function(sigma2,a,rho){  # Get covariance matrix
		# S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		# b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		# P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		# list(P=P,b=b,S=S)
	# }

	# dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # Density of mixture t-distribution
		# x <- matrix(x,,2)	
		# y <- matrix(y,,2)
		# if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# if(!is.null(lc)){
			# S <- S[[lc]]
			# Q <- Q[[lc]]
		# }
		# P <- ifelse(z==1,S,Q)[[1]]  # precision matrix
		# b <- ifelse(z==1,S$b,Q$b)  # determinant
		# d <- x-y
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# out
	# }

	# # Test dmvt2 function
	# # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # test2 <- log(0.5)+
		# # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # plot(test1,test2)
	# # summary(test2-test1)
	# # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # test2 <- 0.5*
		# # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # plot(test1,test2)
	# # summary(test2/test1)
	# # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # Accept/reject proposals for mu
		# mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# num.z1 <- denom.z1 <- num.z0 <- denom.z0 <- 0
		# if(length(idx.1)>0){
			# num.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy.star,z=1,Sigma,Q,lc[x],log=TRUE)))
			# denom.z1 <- sum(sapply(idx.1,function(x)
				# dt2(s[x,],mu.xy,z=1,Sigma,Q,lc[x],log=TRUE)))
		# }
		# if(length(idx.0)>0){
			# num.z0 <- sum(sapply(idx.0,function(x) 
				# dt2(s[x,],mu.xy.star,z=0,Sigma,Q,lc[x],log=TRUE)))
			# denom.z0 <- sum(sapply(idx.0,function(x)
				# dt2(s[x,],mu.xy,z=0,Sigma,Q,lc[x],log=TRUE)))
		# }
		# mh.star <- num.z1+num.z0+U[mu.star[x],]%*%gamma  # numerator of Metropolis ratio
			# #-log(sum(exp(U%*%gamma))))  # integral over S.tilde
			# # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# mh.0 <-	denom.z1+denom.z0+U[mu[x],]%*%gamma   # denominator of Metropolis ratio
			# #-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# }

	# #####################################################################
	# ###  Get starting values from previous model
	# #####################################################################

	# k <- nrow(start$sigma)
	# if(!is.null(k)){ 
		# start$sigma <- start$sigma[k,]
		# start$a <- start$a[k,]
		# start$rho <- start$rho[k,]
		# start$theta <- start$theta[k]+0.1				
		# start$sigma.mu <- start$sigma.mu[k]
		# start$gamma <- start$gamma[k,] 
		# start$h <- start$mu[,k]
		# start$h <- S.tilde[match(start$h,S.tilde[,2]),1]
		# start$z <- start$z[,k]
	# }


	# #####################################################################
	# ###  Setup Variables 
	# #####################################################################
# # browser() 
	# cat("\nSetting up variables....")
	# Ts <- nrow(s)  # number of telemetry locations
	# #* Ty <- length(y)  # number of wet/dry observations
	# Ty <- 0
	# qX <- ncol(X)
	# #* qW <- ncol(W)
	# qW <- 0
	# qU <- ncol(U)
	# v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# #* y1 <- which(y==1)+Ts
	# #* y0 <- which(y==0)+Ts
	# #* y1.sum <- length(y1)
	# #* y0.sum <- length(y0)
	# #* W.cross <- t(W)%*%W  # cross product of W
	# lc <- as.numeric(lc)  # Argos location quality class
	# n.lc <- length(unique(lc))  # number of error classes
	# lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	

	# #####################################################################
	# ### Priors
	# #####################################################################
	
	# cat("\nGetting priors....")
	
	# # Observation model
	# u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	# mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	# sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	# #* r.sigma.alpha <- priors$r.sigma.alpha  # IG prior on sigma.alpha
	# #* q.sigma.alpha <- priors$q.sigma.alpha  # IG prior on sigma.alpha
	# r.theta <- priors$r.theta  # IG prior on theta
	# q.theta <- priors$q.theta  # IG prior on theta
	# sigma.gamma <- priors$sigma.gamma  # variance of normal prior on gamma
	
	# # Temporal haul-out process model 
	# #* mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # Spatial haul-out process model 
	# J <- priors$J  # maximum number of clusters per truncation approximation
	# mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# # Sigma.gamma <- diag(qU)*sigma.gamma^2
  	# # Sigma.gamma.inv <- solve(Sigma.gamma)


	# #####################################################################
	# ### Standardize parameters
	# #####################################################################

	# cat("\nStandardizing variables....")
	
	# # Center and scale s and S.tilde
	# s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.sd
	
	# # Center and scale tuning parameters
	# tune$sigma.mu <- tune$sigma.mu/s.sd
	# tune$mu <- tune$mu/s.sd
	# tune$sigma <- tune$sigma/s.sd

	# # Center and scale priors
	# mu.sigma <- mu.sigma/s.sd  
	# u.sigma <- u.sigma/s.sd

	# # Center and scale starting values
	# sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# sigma <- start$sigma/s.sd  # observation model standard deviation


	# #####################################################################
	# ### Starting values: Appendix A, Steps 1 and 2
	# #####################################################################

	# cat("\nGetting starting values....")
# # browser()
	# # Observation model
	# a <- start$a
	# rho <- start$rho
	# Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	# Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# # Temporal haul-out process model
	# #* alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	# #* Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# #* Sigma.alpha.inv <- solve(Sigma.alpha)
  	# #* A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# z <- start$z  # haul-out status; z=1: dry, z=0: wet

	# # Spatial haul-out process model
	# gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# theta <- start$theta  # DP concentration parameter


  	# #####################################################################
	# ### Create receptacles for output
	# #####################################################################
  
  	# cat("\nCreating receptacles for output....")
	# sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudianl telemetry measurement error
	# a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	# rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
	# beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# #* alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# m.save <- numeric(n.mcmc)  # number of clusters
	# v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# #* sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects
	
    
	# #####################################################################
	# ### MCMC loop: Appendix A, Steps 3-9
	# #####################################################################
	
	# # Track overall MH accpetance rate
	# keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# # Track MH accpetance rate for adaptive tuning
	# keep.tmp <- keep  
	# T.b <- 50  # frequency of adaptive tuning
	# m.save.tmp <- 0  # number of clusters

	# t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# t.mcmc.start <- Sys.time()  # timing MCMC iterations
	
	# cat("\nEntering MCMC Loop....\n")
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# } 
# # browser()

		# if(adapt==TRUE & k%%50==0) {  # Adaptive tuning
			# # browser()
			# keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# tune$sigma.mu <- get.tune(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# tune$gamma <- get.tune(tune$gamma,keep.tmp$gamma,k)
			# tune$mu <- get.tune(tune$mu,keep.tmp$mu,k)
			# tune$sigma <- sapply(1:n.lc,function(x) 
				# get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			# tune$a <- sapply(1:n.lc,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			# tune$rho <- sapply(1:n.lc,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			# keep.tmp <- lapply(keep.tmp,function(x) x*0)
			# m.save.tmp <- 0
	   	# } 	
		
	    
	    # #--------------------------------------------------------------------------
	  	# # Update sigma.mu: Appendix A, Step 4 
	    # #--------------------------------------------------------------------------
# # browser()
		# # Lognormal prior
	    # sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		# if(sigma.mu.star>0){
			# Q.star <- sapply(1:n.lc,function(x) 
				# get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
			# idx <- which(z==0)
		    # mh.star.sigma.mu <-	sum(sapply(idx,function(x) 
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    # mh.0.sigma.mu <- sum(sapply(idx,function(x)
		    	# dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	# dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# sigma.mu <- sigma.mu.star
				# Q <- Q.star
				# keep$sigma.mu <- keep$sigma.mu+1
				# keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	# } 
		# }

	
		# #--------------------------------------------------------------------------
		# # Update observation model parameters: Appendix A, Step 5
		# #--------------------------------------------------------------------------
# # browser()
		# sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		# a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		# rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		# for(i in 1:n.lc){  # loop to iterate over error classes: Appendix A, Step 5(b)

			# idx <- lc.list[[i]]  # index of locations in error class i
			# z1 <- idx[which(z[idx]==1)]
			# z0 <- idx[which(z[idx]==0)]

			# ###
			# ### Sample sigma: Appendix A, Step 5(a.i)
			# ###

			# if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				# Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# mh.star.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# mh.0.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# sigma[i] <- sigma.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$sigma[i] <- keep$sigma[i]+1
					# keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				# }
			# }

			# ###
			# ### Sample a: Appendix A, Step 5(a.ii)
			# ###
			
			# if(a.star[i]>0 & a.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# mh.star.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# mh.0.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.a-mh.0.a)>runif(1)){
					# a[i] <- a.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$a[i] <- keep$a[i]+1
					# keep.tmp$a[i] <- keep.tmp$a[i]+1
				# }
			# }

			# ###
			# ### Sample rho: Appendix A, Step 5(a.iii)
			# ###
			
			# if(rho.star[i]>0 & rho.star[i]<1){
				# Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# mh.star.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# mh.0.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# rho[i] <- rho.star[i]
					# Sigma[[i]] <- Sigma.star
					# Q[[i]] <- Q.star
					# keep$rho[i] <- keep$rho[i]+1
					# keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				# }
			# }
		# }

	
		# #--------------------------------------------------------------------------
		# # Update temporal haul-out process model parameters: Appendix A, Step 6
		# #--------------------------------------------------------------------------
# # browser()
	 	# t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# ###	
		# ### Update beta: Appendix A, Step 6(a)
		# ###
		
		# #* b <- crossprod(X,(v-W%*%alpha))
		# b <- crossprod(X,v)
	  	# beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		# ###
		# ### Update alpha: Appendix A, Step 6(b)
		# ###
		
		# #* b <- crossprod(W,(v-X%*%beta))
		# #* alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# ###
		# ### Update v(t_y): Appendix A, Step 6(c)
		# ###

		# linpred <- X%*%beta  #* +W%*%alpha  # update linear predictor 
	  	# #* v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	# #* v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# ###
		# ### Update v(t_s): Appendix A, Step 6(d)
		# ###
			
		# z1 <- which(z==1)
		# z0 <- which(z==0)
		# v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	# v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	# t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of aux. variable update
	  	
	  	# ###
	  	# ### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	# ###
	  	
		# #* p <- pnorm(linpred[1:Ts,])

	    # ###
	    # ### Update z: Appendix A, Step 6(f)
	    # ###
		
		# #* p1 <- p*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		# #* p2 <- (1-p)*sapply(1:Ts,function(x) 
			# # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		# #* p <- exp(log(p1)-log(p1+p2))	
		# #* z <- rbinom(Ts,1,p)

		# ###
		# ### Update sigma2.alpha: Appendix A, Step 6(g)
		# ###
		
		# #* r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.sigma.alpha)
		# #* q.tmp <- qW/2+q.sigma.alpha
		# #* sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# #* diag(Sigma.alpha) <- sigma2.alpha
		# #* Sigma.alpha.inv <- solve(Sigma.alpha)
		# #* A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		
					
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 7
		# #--------------------------------------------------------------------------
# # browser()
		# ###
		# ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# ###

		# n <- table(h)  
		# m <- length(n)  # number of clusters
		# mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		# ###
		# ### Update the stick-breaking process: Appendix A, Step 7(b)
		# ###
		
			# # Create index set I: Appendix A, Step 6(b.i)
			# I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # Update eta: Appendix A, Step 6(b.ii)
			# n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # Update pi: Appendix A, Step 6(b.ii)
		    # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # Update theta: Appendix A, Step 6(b.iii)
			# theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  

		# ###
		# ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# ###

		# p <- exp(U[-mu,]%*%gamma)
		# mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu
				
		# ###
		# ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# ###
# # browser()			
		# mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# keep$mu <- keep$mu+sum(mh)
		# keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# ###
	    # ### Update gamma: Appendix A, Step 7(e)
	    # ###
# # browser()	    
		# # Integral over S.tilde, occupied mu only
		# gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,sigma.gamma,log=TRUE))+
			# sum(n*c(U[mu,]%*%gamma.star)-n*log(sum(exp(U%*%gamma.star))))
			# # sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
		# mh.0.gamma <- sum(dnorm(gamma,mu.gamma,sigma.gamma,log=TRUE))+
			# sum(n*c(U[mu,]%*%gamma)-n*log(sum(exp(U%*%gamma))))
			# # sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # gamma <- gamma.star
	        # # gamma.int <- gamma.int.star
	        # keep$gamma <- keep$gamma+1
	        # keep.tmp$gamma <- keep.tmp$gamma+1
    	# } 
			
		# ###
	    # ### Update h(t_s): Appendix A, Step 7(f)
	    # ###

		# mu <- c(mu[I],mu.tmp)
		# h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# #--------------------------------------------------------------------------
		# # Save samples: Appendix A, Step 8
		# #--------------------------------------------------------------------------
						
		# sigma.save[k,] <- sigma
		# a.save[k,] <- a
		# rho.save[k,] <- rho
		# mu.save[,k] <- S.tilde[h,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# #* sigma.alpha.save[k] <- sigma2.alpha
		# #* alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# gamma.save[k,] <- gamma
		# v.save[,k] <- v
		# z.save[,k] <- z
		# m.save[k] <- m
		# m.save.tmp <- m.save.tmp+m
	# }  # end of MCMC loop

  	# t.mcmc.end <- Sys.time()
  	
  	# tune$sigma.mu <- tune$sigma.mu*s.sd
  	# tune$sigma <- tune$sigma*s.sd
	# tune$mu <- tune$mu*s.sd
  	# sigma.save <- sigma.save*s.sd
  	# sigma.mu.save <- sigma.mu.save*s.sd
	# #* sigma.alpha.save <- sqrt(sigma.alpha.save)
	
  	# #####################################################################
	# ### Write output
	# #####################################################################
	  
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu <- keep$mu/sum(m.save)
	# keep$gamma <- keep$gamma/n.mcmc
	# keep$sigma <- keep$sigma/n.mcmc
	# keep$a <- keep$a/n.mcmc
	# keep$rho <- keep$rho/n.mcmc

	# cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	# cat("\na acceptance rate:",round(keep$a,2))
	# cat("\nrho acceptance rate:",round(keep$rho,2)) 

	# cat(paste("\n\nEnd time:",Sys.time()))
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# cat(paste("\nTime per MCMC iteration:",
		# round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	# list(beta=beta.save,gamma=gamma.save,  #* alpha=alpha.save,
		# mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# sigma.mu=sigma.mu.save,  #* sigma.alpha=sigma.alpha.save,
		# sigma=sigma.save,a=a.save,rho=rho.save,
		# keep=keep,tune=tune,n.mcmc=n.mcmc)
# }









# # 
# haulouts.mcmc <- function(s,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# ###
 	# ### Brian M. Brost (08 DEC 2015)
 	# ### Cleaned up version of haulouts.3.mcmc.R
 	# ###
 	
 	# ###
 	# ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# ### the location of haul-out sites; S.tilde=support of haul-out sites 
 	# ###

	# t.start <- Sys.time()
	# cat(paste("Start time:",t.start,"\n"))

	# #####################################################################
	# ### Libraries and Subroutines
	# #####################################################################
  
	# # library(MCMCpack)  # for Dirichlet distribution functions
	# # library(data.table)  # for tabulating and summing
	# # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# # library(doParallel)  # for parallel processing
	# # library(foreach)  # for parallel processing  
	# # library(mvtnorm)  # for multivariate normal density
	# # library(msm)  # for truncated normal density
  
  	# truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# }

	# adapt <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		# exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# }
	
	# dt2 <- function(x,y,z,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # browser()		
		# x <- matrix(x,,2)	
		# y <- matrix(y,,2)
		# if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# P <- ifelse(z==1,Sigma,Q)[[1]]  # precision matrix
		# b <- ifelse(z==1,Sigma$b,Q$b)  # determinant
		# d <- x-y
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# out
	# }

	# get.Sigma <- function(sigma2,a,rho){  # get var-cov matrix, determinant, and inverse
		# n.lc <- length(sigma2)
		# Sigma <- lapply(1:n.lc,function(x) sigma2[x]*  # variance-covariance matrix			
			# matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))  
		# det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		# P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			# matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		# list(Sigma=Sigma,P=P,det=det)  # return list of lists
	# }

	
	
	
	# dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # Calculate density of mixture t distribution
		# x <- matrix(x,,2,byrow=FALSE)
		# if(nrow(x)==0) out <- 0
		# if(nrow(x)!=0){
			# if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# lc.idx <- sort(unique(lc))
			# n <- length(lc.idx)
			# lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# P <- Sigma$P[lc.idx]  # precision matrix
			# b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))  # determinant
			# d <- x-y  
			# out <- numeric(nrow(d))
			# for(i in 1:n){  # calculate kernel for each Sigma
				# idx <- lc.list[[i]]
				# d.tmp <- d[idx,]
				# out[idx] <- tkern(d.tmp,P[[i]],nu)+tkern(d.tmp,K%*%P[[i]]%*%t(K),nu)
			# }	
			# if(log) out <- log(0.5)+log(out)+log(b[lc])  # log density
					# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			# if(!log){ out <- 0.5*out*b[lc]  # density
					# # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
			# } 
		# }
		# out
	# }

	# # Test dmvt2 function
	# # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # test2 <- log(0.5)+
		# # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # plot(test1,test2)
	# # summary(test2-test1)
	# # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # test2 <- 0.5*
		# # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # plot(test1,test2)
	# # summary(test2/test1)
	# # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # browser()
		# # Accept/reject proposals for mu
		# mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# mh.star <- sum(  # numerator of Metropolis-Hastings ratio
			# dmvt2(s[idx.1,],mu.xy.star,lc[idx.1],Sigma,log=TRUE),
			# dmvt2(s[idx.0,],mu.xy.star,lc[idx.0],Q,log=TRUE))+
			# U[mu[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# mh.0 <-	sum(  # denominator of Metropolis-Hastings ratio
			# dmvt2(s[idx.1,],mu.xy,lc[idx.1],Sigma,log=TRUE),
			# dmvt2(s[idx.0,],mu.xy,lc[idx.0],Q,log=TRUE))+
			# U[mu.star[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# }

	# get.h <- function(x,y,z,lc,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2)){
		# # browser()
		# # For sampling h
		# if(z==1) P <- Sigma$P[[lc]]  # precision matrix when z=1 
		# if(z==0) P <- Q$P[[lc]]  # precision matrix when z=0
		# d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# log(out)  # log of mixture kernel
	# }


	# #####################################################################
	# ###  Setup Variables 
	# #####################################################################
# # browser() 
	# cat("\nSetting up variables....")

	# Ts <- nrow(s)  # number of telemetry locations
	# # Ty <- length(y)  # number of wet/dry observations
# Ty <- 0
	# X <- as.matrix(X)
	# U <- as.matrix(U)
	# qX <- ncol(X)
	# # qW <- ncol(W)
# qW <- 0
	# # qU <- dim(U)[3]+1	# number of raster layers plus an intercept
# qU <- ncol(U)
	# v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# # y1 <- which(y==1)+Ts
	# # y0 <- which(y==0)+Ts
	# # y1.sum <- length(y1)
	# # y0.sum <- length(y0)
	# # W.cross <- t(W)%*%W  # cross product of W
	# # idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	# # S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# # information in S.tilde; note that mu below references row idx in S.tilde
	# h <- match(start$h,S.tilde[,2])  # Note: h corresponds to row idx of S.tilde 
		# # for computational efficiency, and not idx of mu as in Appendix A
	# # U <- cbind(1,values(U)[idx])  # convert raster to design matrix

	# lc <- as.numeric(priors$lc)  # Argos location quality class
	# n.lc <- length(unique(lc))  # number of error classes
	# lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	
	# #####################################################################
	# ### Standardize parameters
	# #####################################################################

	# cat("\nStandardizing variables....")
	
	# # Center and scale s and S.tilde
	# s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# ncol=2,byrow=TRUE))/s.sd
	
	# # Center and scale tuning parameters
	# tune$sigma.mu <- tune$sigma.mu/s.sd
	# tune$mu <- tune$mu/s.sd
	# tune$sigma <- tune$sigma/s.sd

	# # Center and scale priors
	# priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	# priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu
	# priors$mu.sigma <- priors$mu.sigma/s.sd  # variance on lognormal prior for sigma.mu

	# # Center and scale starting values
	# sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# sigma <- start$sigma/s.sd  # observation model standard deviation


	# #####################################################################
	# ### Priors
	# #####################################################################
	
	# cat("\nGetting priors....")
	
	# # Observation model
	# u.sigma <- priors$u.sigma/s.sd
	
	# # Temporal haul-out process model 
	# # mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# Sigma.beta.inv <- solve(Sigma.beta)
	# A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # Spatial haul-out process model 
	# J <- priors$J  # maximum number of clusters per truncation approximation
	# mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	# Sigma.gamma.inv <- solve(Sigma.gamma)


	# #####################################################################
	# ### Appendix A, Step 1: starting values 
	# #####################################################################

	# cat("\nGetting starting values....")

	# # Observation model
	# a <- start$a
	# rho <- start$rho
	# Sigma <- get.Sigma(sigma^2,a,rho)
	# Q <- get.Sigma(sigma^2+sigma.mu^2,a,rho)

	# # Temporal haul-out process model
	# beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	# # alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	# # Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# # Sigma.alpha.inv <- solve(Sigma.alpha)
  	# # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # linpred <- X%*%beta+W%*%alpha
# linpred <- X%*%beta

	# # Spatial haul-out process model
	# gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# theta <- start$theta  # DP concentration parameter


  	# #####################################################################
	# ### Create receptacles for output
	# #####################################################################
  
  	# cat("\nCreating receptacles for output....")
	# beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# # alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# m.save <- numeric(n.mcmc)  # number of clusters
	# v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

	# sigma.save <- matrix(0,n.mcmc,n.lc)
	# a.save <- matrix(0,n.mcmc,n.lc)
	# rho.save <- matrix(0,n.mcmc,n.lc)

    
	# #####################################################################
	# ### Appendix A, Steps 3-8: MCMC loop 
	# #####################################################################

	# keep <- list(mu=0,sigma.mu=0,gamma=0)  # number of proposals accepted for Metropolis updates

# keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))
	# keep.tmp <- list(mu=0,sigma.mu=0,gamma=0)  # for adaptive tuning
	# m.save.tmp <- 0
	# t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# t.mcmc.start <- Sys.time()  # timing MCMC iterations
	# T.b <- 50  # frequency of adaptive tuning
	
	# cat("\nEntering MCMC Loop....\n")
	# for (k in 1:n.mcmc) {
    	# if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# cat(k,"");flush.console()	
    		# plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# } 

    	# if(k%%50==0) {  # Adaptive tuning
			# # browser()
			# keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# tune$sigma.mu <- adapt(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# tune$gamma <- adapt(tune$gamma,keep.tmp$gamma,k)
			# tune$mu <- adapt(tune$mu,keep.tmp$mu,k)
			# keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	# } 
	
		# #--------------------------------------------------------------------------
		# # Appendix A, Step 4: update observation model parameters (Sigma)
		# #--------------------------------------------------------------------------

# # browser()
# # priors$sigma <- priors$sigma/s.scale

		# sigma.star <- rnorm(n.lc,sigma,tune$sigma) #Proposals for sigma
		# a.star <- rnorm(n.lc,a,tune$a) #Proposals for a
		# rho.star <- rnorm(n.lc,rho,tune$rho) #Proposals for rho

	# dmvt2.vec <- function(d,P,b,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # Calculate density of mixture t distribution
		# b <- b^(-0.5)  # determinant
		# out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# out <- log(0.5)+log(out)+log(b)  # log density
					# # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# out
	# }
	
		# for(i in 1:n.lc){ #Loop to iterate over error classes: Appendix A, step 2(f)
	# # i <- 2
			# idx <- lc.list[[i]] #Index of locations in error class i
			# d.tmp <- s[idx,]-S.tilde[h[idx],3:4]
			# K <- matrix(c(-1,0,0,1),2)
			# z1 <- which(z[idx]==1)
			# z0 <- which(z[idx]==0)
# # browser()	
			# ### Sample sigma: Appendix A, step 2(b)

			# if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])

				# mh.star.sigma <- 
				# sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# mh.0.sigma <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))

		# # mh.star.sigma <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		# # mh.0.sigma <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))

				# if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# sigma[i] <- sigma.star[i]
					# keep$sigma[i] <- keep$sigma[i]+1
					# Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# Sigma$P[[i]] <- Sigma.star$P[[1]]
					# Sigma$det[[i]] <- Sigma.star$det[[1]]
					# Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# Q$P[[i]] <- Q.star$P[[1]]
					# Q$det[[i]] <- Q.star$det[[1]]
				# }
			# }

			# ### Sample a: Appendix A, step 2(c)

			# if(a.star[i]>0 & a.star[i]<1){
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])

				# mh.star.a <- 
				# sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# mh.0.a <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))
				# # mh.star.a <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
			  	# # mh.0.a <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
				# if(exp(mh.star.a-mh.0.a)>runif(1)){
					# a[i] <- a.star[i]
					# keep$a[i] <- keep$a[i]+1
					# Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# Sigma$P[[i]] <- Sigma.star$P[[1]]
					# Sigma$det[[i]] <- Sigma.star$det[[1]]
					# Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# Q$P[[i]] <- Q.star$P[[1]]
					# Q$det[[i]] <- Q.star$det[[1]]
				# }
			# }

			# ### Sample rho: Appendix A, step 2(d)

			# if(rho.star[i]>0 & rho.star[i]<1){
				# Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])

				# mh.star.rho <- 
				# sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# mh.0.rho <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))

				# # mh.star.rho <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
			  	# # mh.0.rho <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
				# if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# rho[i] <- rho.star[i]
					# keep$rho[i] <- keep$rho[i]+1
					# Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# Sigma$P[[i]] <- Sigma.star$P[[1]]
					# Sigma$det[[i]] <- Sigma.star$det[[1]]
					# Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# Q$P[[i]] <- Q.star$P[[1]]
					# Q$det[[i]] <- Q.star$det[[1]]
				# }
			# }
		# }


	
		# #--------------------------------------------------------------------------
		# # Appendix A, Step 4: update temporal haul-out process model parameters 
		# #--------------------------------------------------------------------------
# # browser()
	 	# t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# ###
		# ### Appendix A, Step 4(a): update v(t_y)
		# ###
		
	  	# # v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# # v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# ###
		# ### Appendix A, Step 4(b): update v(t_s)
		# ###
			
		# z1 <- which(z==1)
		# z0 <- which(z==0)
	  	# v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	
		# ###
		# ### Appendix A, Step 4(c): update alpha
		# ###
		
		# # b <- crossprod(W,(v-X%*%beta))
		# # alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# ###	
		# ### Appendix A, Step 4(d): update beta
		# ###
		
		# # b <- crossprod(X,(v-W%*%alpha))
# # browser()
# b <- crossprod(X,v)
	  	# beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	# t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of auxilliary variable update
	  	
	  	# ###
	  	# ### Appendix A, Step 4(e): calculate Prob(z(t_s)==1)
	  	# ###
	  	
	  	# # linpred <- X%*%beta+W%*%alpha  # update linear predictor 
		# # p <- pnorm(linpred[1:Ts,])

	    # ###
	    # ### Appendix A, Step 4(f): update z
	    # ###
	    
		# # p1 <- p*dmvt2(s,S.tilde[h,3:4],lc,Sigma,log=FALSE)
		# # p2 <- (1-p)*dmvt2(s,S.tilde[h,3:4],lc,Q,log=FALSE)
		# # p <- exp(log(p1)-log(p1+p2))	
		# # z <- rbinom(Ts,1,p)

		# ###
		# ### Appendix A, Step 4(g): update sigma2.alpha
		# ###
		
		# # r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		# # q.tmp <- qW/2+priors$q.sigma.alpha
		# # sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# # diag(Sigma.alpha) <- sigma2.alpha
		# # Sigma.alpha.inv <- solve(Sigma.alpha)
		# # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		

	    # #--------------------------------------------------------------------------
	  	# # Appendix A, Step 5: update sigma.mu
	    # #--------------------------------------------------------------------------
# # browser()

		# # test <- rnorm(1000,log(mu.sigma),tau)
		# # hist(s.sd*exp(test))

		# # Lognormal prior
	    # sigma.mu.star <-  exp(rnorm(1,log(sigma.mu),tune$sigma.mu))
		# Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,a,rho)
		# idx <- which(z==0)
	    # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))					+dnorm(log(sigma.mu.star),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))						    +dnorm(log(sigma.mu),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	# sigma.mu <- sigma.mu.star
			# Q <- Q.star
			# keep$sigma.mu <- keep$sigma.mu+1
			# keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1

    	# } 

		# # Uniform prior	    
	    # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)				   
	    # # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# # # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,n.lc,Mix)
			# # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,a,rho)
			# # idx <- which(z==0)
		    # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		    # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
		    # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# # sigma.mu <- sigma.mu.star
				# # Q <- Q.star
				# # keep$sigma.mu <- keep$sigma.mu+1
	    	# # } 
	    # # }
		
				
		# #--------------------------------------------------------------------------
		# # Update spatial haul-out process model parameters: Appendix A, Step 6
		# #--------------------------------------------------------------------------
# # browser()
		# ###
		# ### Appendix A, Step 6(a): tabulate cluster membership 
		# ###

		# n <- table(h)  
		# m <- length(n)  # number of clusters
		# mu <- as.numeric(names(n))  # idx of occupied clusters

		# ###
		# ### Appendix A, Step 6(d): update 'unoccupied' mu
		# ###
# # browser()
		# p <- exp(U[-mu,]%*%gamma)
		# mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		# ###
		# ### Appendix A, Step 6(c): update 'occupied' mu		   
		# ###
# # browser()			
		# mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# mh <- sapply(dup.idx,function(x)  # accepted proposals
			# get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# keep$mu <- keep$mu+sum(mh)
		# keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# keep.idx <- dup.idx[mh]
		# mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# ###
	    # ### Appendix A, Step 6(f): update gamma
	    # ###

		# # Integral over S.tilde, occupied mu only
		# gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
			# mh.0.gamma <- sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # gamma <- gamma.star
	        # # gamma.int <- gamma.int.star
	        # keep$gamma <- keep$gamma+1
	        # keep.tmp$gamma <- keep.tmp$gamma+1
    	# } 

		# ###
		# ### Appendix A, Step 6(b): update the stick-breaking process 
		# ###
		
			# # Appendix A, Step 6(b(i)): create index set I
			# I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # Appendix A, Step 6(b(ii)): update eta
			# n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # Appendix A, Step 6(b(ii)): update pi
		    # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # Appendix A, Step 6(b(iii)): update theta
			# theta <- rgamma(1,priors$r.theta+J-1,priors$q.theta-sum(log(1-eta[-J])))  

			
		# ###
	    # ### Appendix A, Step 6(e): update h(t_s)
	    # ###

		# mu <- c(mu[I],mu.tmp)
		# h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# exp(log(pie)+get.h(s[x,],S.tilde[mu,3:4],z[x],lc[x],Sigma,Q))))

		# ###
		# ###  Appendix A, Step 7: save samples 		   
		# ###
		
		# mu.save[,k] <- S.tilde[h,2]
		# theta.save[k] <- theta    
		# sigma.mu.save[k] <- sigma.mu
		# # sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		# # alpha.save[k,] <- alpha
		# beta.save[k,] <- beta
		# gamma.save[k,] <- gamma
		# v.save[,k] <- v
		# z.save[,k] <- z
		# m.save[k] <- m
		# m.save.tmp <- m.save.tmp+m
# sigma.save[k,] <- sigma*s.sd
# a.save[k,] <- a
# rho.save[k,] <- rho
	# }
  	
  	# tune$sigma.mu <- tune$sigma.mu*s.sd
	# tune$mu <- tune$mu*s.sd
  	# sigma.mu.save <- sigma.mu.save*s.sd
  	# t.mcmc.end <- Sys.time()

	# #####################################################################
	# ### Write output
	# #####################################################################
	  
	# keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# keep$mu <- keep$mu/sum(m.save)
	# keep$gamma <- keep$gamma/n.mcmc
	# keep$sigma <- keep$sigma/n.mcmc
	# keep$a <- keep$a/n.mcmc
	# keep$rho <- keep$rho/n.mcmc

	# cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 

	# cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	# cat(paste("\na acceptance rate:",round(keep$a,2))) 
	# cat(paste("\nrho acceptance rate:",round(keep$rho,2))) 

	# cat(paste("\n\nEnd time:",Sys.time()))
	# cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# cat(paste("\nTime per MCMC iteration:",
		# round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	# list(beta=beta.save,gamma=gamma.save,  # alpha=alpha.save,
		# mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# sigma.mu=sigma.mu.save,  # sigma.alpha=sigma.alpha.save,
		# keep=keep,tune=tune,n.mcmc=n.mcmc,
		# sigma=sigma.save,a=a.save,rho=rho.save)
# }



# # haulouts.mcmc <- function(s,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# # ###
 	# # ### Brian M. Brost (08 DEC 2015)
 	# # ### See haulouts.sim.R to simulate data according to this model specification,
 	# # ### and haulouts.script.R for implementations with harbor seal data
 	# # ###
 	
 	# # ###
 	# # ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# # ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# # ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# # ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# # ### the location of haul-out sites; S.tilde=support of haul-out sites 
 	# # ###

	# # t.start <- Sys.time()
	# # cat(paste("Start time:",t.start,"\n"))

	# # #####################################################################
	# # ### Libraries and Subroutines
	# # #####################################################################
  
	# # # library(MCMCpack)  # for Dirichlet distribution functions
	# # # library(data.table)  # for tabulating and summing
	# # # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# # # library(doParallel)  # for parallel processing
	# # # library(foreach)  # for parallel processing  
	# # # library(mvtnorm)  # for multivariate normal density
	# # # library(msm)  # for truncated normal density
  
  	# # truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# # flow <- pnorm(low,mu,sqrt(sig2)) 
		# # fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# # u <- runif(nsamp) 
		# # tmp <- flow+u*(fhigh-flow)
		# # x <- qnorm(tmp,mu,sqrt(sig2))
		# # x
	# # }
	
	# # tkern <- function(d,P,nu){  # kernel of t-distribution
		# # (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# # }

	# # get.Sigma <- function(sigma2,a,rho){  # get var-cov matrix, determinant, and inverse
		# # n.lc <- length(sigma2)
		# # Sigma <- lapply(1:n.lc,function(x) sigma2[x]*  # variance-covariance matrix			
			# # matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))  
		# # det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		# # P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			# # matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		# # list(Sigma=Sigma,P=P,det=det)  # return list of lists
	# # }

	# # dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # # Calculate density of mixture t distribution
		# # x <- matrix(x,,2,byrow=FALSE)
		# # if(nrow(x)==0) out <- 0
		# # if(nrow(x)!=0){
			# # if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# # lc.idx <- sort(unique(lc))
			# # n <- length(lc.idx)
			# # lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# # P <- Sigma$P[lc.idx]  # precision matrix
			# # b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))  # determinant
			# # d <- x-y  
			# # out <- numeric(nrow(d))
			# # for(i in 1:n){  # calculate kernel for each Sigma
				# # idx <- lc.list[[i]]
				# # d.tmp <- d[idx,]
				# # out[idx] <- tkern(d.tmp,P[[i]],nu)+tkern(d.tmp,K%*%P[[i]]%*%t(K),nu)
				# # # out[idx] <- (1+1/nu*(rowSums((d.tmp%*%P[[i]])*d.tmp)))^-((nu+2)/2) +
					# # # (1+1/nu*(rowSums((d.tmp%*%(K%*%P[[i]]%*%t(K)))*d.tmp)))^-((nu+2)/2)
			# # }	
			# # if(log) out <- log(0.5)+log(out)+log(b[lc])  # log density
					# # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			# # if(!log){ out <- 0.5*out*b[lc]  # density
					# # # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
			# # } 
		# # }
		# # out
	# # }

	# # # Test dmvt2 function
	# # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # # test2 <- log(0.5)+
		# # # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # plot(test1,test2)
	# # # summary(test2-test1)
	# # # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # # test2 <- 0.5*
		# # # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # plot(test1,test2)
	# # # summary(test2/test1)
	# # # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# # get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # # browser()
		# # # Accept/reject proposals for mu
		# # mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# # mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# # idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# # idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# # mh.star <- sum(  # numerator of Metropolis-Hastings ratio
			# # dmvt2(s[idx.1,],mu.xy.star,lc[idx.1],Sigma,log=TRUE),
			# # dmvt2(s[idx.0,],mu.xy.star,lc[idx.0],Q,log=TRUE))+
			# # U[mu[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# # mh.0 <-	sum(  # denominator of Metropolis-Hastings ratio
			# # dmvt2(s[idx.1,],mu.xy,lc[idx.1],Sigma,log=TRUE),
			# # dmvt2(s[idx.0,],mu.xy,lc[idx.0],Q,log=TRUE))+
			# # U[mu.star[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# # exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# # }

	# # get.h <- function(x,y,z,lc,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2)){
		# # # browser()
		# # # For sampling h
		# # if(z==1) P <- Sigma$P[[lc]]  # precision matrix when z=1 
		# # if(z==0) P <- Q$P[[lc]]  # precision matrix when z=0
		# # d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# # # out <- (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2) +
			# # # (1+1/nu*(rowSums((d%*%(K%*%P%*%t(K)))*d)))^-((nu+2)/2)
		# # log(out)  # log of mixture kernel
	# # }


	# # #####################################################################
	# # ###  Setup Variables 
	# # #####################################################################
# # # browser() 
	# # cat("\nSetting up variables....")

	# # Ts <- nrow(s)  # number of telemetry locations
	# # # Ty <- length(y)  # number of wet/dry observations
# # Ty <- 0
	# # X <- as.matrix(X)
	# # U <- as.matrix(U)
	# # qX <- ncol(X)
	# # # qW <- ncol(W)
# # qW <- 0
	# # # qU <- dim(U)[3]+1	# number of raster layers plus an intercept
# # qU <- ncol(U)
	# # v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# # # y1 <- which(y==1)+Ts
	# # # y0 <- which(y==0)+Ts
	# # # y1.sum <- length(y1)
	# # # y0.sum <- length(y0)
	# # # W.cross <- t(W)%*%W  # cross product of W
	# # # idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	# # # S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# # # information in S.tilde; note that mu below references row idx in S.tilde
	# # h <- match(start$h,S.tilde[,2])  # Note: h corresponds to row idx of S.tilde 
		# # # for computational efficiency, and not idx of mu as in Appendix A
	# # # U <- cbind(1,values(U)[idx])  # convert raster to design matrix
	
	
	# # #####################################################################
	# # ### Standardize parameters
	# # #####################################################################

	# # cat("\nStandardizing variables....")
	
	# # # Center and scale s and S.tilde
	# # s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# # s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# # s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# # S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# # ncol=2,byrow=TRUE))/s.sd
	
	# # # Center and scale tuning parameters
	# # tune$sigma.mu <- tune$sigma.mu/s.sd
	# # tune$mu <- tune$mu/s.sd

	# # # Center and scale priors
	# # priors$sigma <- priors$sigma/s.sd  # observation model standard deviation
	# # priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	# # priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu

	# # # Center and scale starting values
	# # sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter


	# # #####################################################################
	# # ### Priors
	# # #####################################################################
	
	# # cat("\nGetting priors....")
	
	# # # Observation model
	# # lc <- as.numeric(priors$lc)  # Argos location quality class
	# # sigma2 <- priors$sigma^2  # standardized variance component
	# # Sigma <- get.Sigma(sigma2,priors$a,priors$rho)
	# # Q <- get.Sigma(sigma2+sigma.mu^2,priors$a,priors$rho)
# # # browser()
	# # # Temporal haul-out process model 
	# # mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# # mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# # Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# # Sigma.beta.inv <- solve(Sigma.beta)
	# # A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # # Spatial haul-out process model 
	# # J <- priors$J  # maximum number of clusters per truncation approximation
	# # mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# # Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	# # Sigma.gamma.inv <- solve(Sigma.gamma)


	# # #####################################################################
	# # ### Appendix A, Step 1: starting values 
	# # #####################################################################

	# # cat("\nGetting starting values....")

	# # # Temporal haul-out process model
	# # beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	# # # alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	# # # Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# # # Sigma.alpha.inv <- solve(Sigma.alpha)
  	# # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# # z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # # linpred <- X%*%beta+W%*%alpha
# # linpred <- X%*%beta
	# # # Spatial haul-out process model
	# # gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# # theta <- start$theta  # DP concentration parameter


  	# # #####################################################################
	# # ### Create receptacles for output
	# # #####################################################################
  
  	# # cat("\nCreating receptacles for output....")
	# # beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# # alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# # gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# # mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# # theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# # m.save <- numeric(n.mcmc)  # number of clusters
	# # v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# # sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

    
	# # #####################################################################
	# # ### Appendix A, Steps 3-8: MCMC loop 
	# # #####################################################################

	# # keep <- list(mu=0,sigma.mu=0,gamma=0)  # number of proposals accepted for Metropolis updates
	# # t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# # t.mcmc.start <- Sys.time()  # timing MCMC iterations

  	# # cat("\nEntering MCMC Loop....\n")
	# # for (k in 1:n.mcmc) {
    	# # if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# # cat(k,"");flush.console()	
    		# # plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# # } 
	
	
		# # #--------------------------------------------------------------------------
		# # # Appendix A, Step 4: update temporal haul-out process model parameters 
		# # #--------------------------------------------------------------------------
# # # browser()
	 	# # t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# # ###
		# # ### Appendix A, Step 4(a): update v(t_y)
		# # ###
		
	  	# # # v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# # # v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# # ###
		# # ### Appendix A, Step 4(b): update v(t_s)
		# # ###
			
		# # z1 <- which(z==1)
		# # z0 <- which(z==0)
	  	# # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	
		# # ###
		# # ### Appendix A, Step 4(c): update alpha
		# # ###
		
		# # # b <- crossprod(W,(v-X%*%beta))
		# # # alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# # ###	
		# # ### Appendix A, Step 4(d): update beta
		# # ###
		
		# # # b <- crossprod(X,(v-W%*%alpha))
# # b <- crossprod(X,v)
	  	# # beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	# # t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of auxilliary variable update
	  	
	  	# # ###
	  	# # ### Appendix A, Step 4(e): calculate Prob(z(t_s)==1)
	  	# # ###
	  	
	  	# # # linpred <- X%*%beta+W%*%alpha  # update linear predictor 
		# # # p <- pnorm(linpred[1:Ts,])

	    # # ###
	    # # ### Appendix A, Step 4(f): update z
	    # # ###
	    
		# # # p1 <- p*dmvt2(s,S.tilde[h,3:4],lc,Sigma,log=FALSE)
		# # # p2 <- (1-p)*dmvt2(s,S.tilde[h,3:4],lc,Q,log=FALSE)
		# # # p <- exp(log(p1)-log(p1+p2))	
		# # # z <- rbinom(Ts,1,p)

		# # ###
		# # ### Appendix A, Step 4(g): update sigma2.alpha
		# # ###
		
		# # # r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		# # # q.tmp <- qW/2+priors$q.sigma.alpha
		# # # sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# # # diag(Sigma.alpha) <- sigma2.alpha
		# # # Sigma.alpha.inv <- solve(Sigma.alpha)
		# # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		

	    # # #--------------------------------------------------------------------------
	  	# # # Appendix A, Step 5: update sigma.mu
	    # # #--------------------------------------------------------------------------
# # # browser()

# # priors$mu.sigma <- 6500/s.sd  # variance on lognormal prior for sigma.mu
# # priors$tau <- 0.1

# # sigma.mu.star <-  exp(rnorm(1,log(sigma.mu),tune$sigma.mu))
		# # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,a,rho)
		# # idx <- which(z==0)
	    # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))					+dnorm(log(sigma.mu.star),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))						    +dnorm(log(sigma.mu),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	# # sigma.mu <- sigma.mu.star
			# # Q <- Q.star
			# # keep$sigma.mu <- keep$sigma.mu+1
			# # # keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1

  	# # } 



	    # # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    # # # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# # # # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,n.lc,Mix)
			# # # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,a,rho)
			# # # idx <- which(z==0)
		    # # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		    # # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
		    # # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# # # sigma.mu <- sigma.mu.star
				# # # Q <- Q.star
				# # # keep$sigma.mu <- keep$sigma.mu+1
	    	# # # } 
	    # # # }
		
				
		# # #--------------------------------------------------------------------------
		# # # Update spatial haul-out process model parameters: Appendix A, Step 6
		# # #--------------------------------------------------------------------------
# # # browser()
		# # ###
		# # ### Appendix A, Step 6(a): tabulate cluster membership 
		# # ###

		# # n <- table(h)  
		# # m <- length(n)  # number of clusters
		# # mu <- as.numeric(names(n))  # idx of occupied clusters

		# # ###
		# # ### Appendix A, Step 6(d): update 'unoccupied' mu
		# # ###

		# # p <- exp(U[-mu,]%*%gamma)
		# # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		# # ###
		# # ### Appendix A, Step 6(c): update 'occupied' mu		   
		# # ###
# # # browser()			
		# # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*	dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# # keep$mu <- keep$mu+sum(mh)
		# # keep.idx <- dup.idx[mh]
		# # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # ###
	    # # ### Appendix A, Step 6(f): update gamma
	    # # ###

		# # # Integral over S.tilde, occupied mu only
		# # gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# # mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
			# # mh.0.gamma <- sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# # if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # # gamma <- gamma.star
	        # # # gamma.int <- gamma.int.star
	        # # keep$gamma <- keep$gamma+1
    	# # } 

		# # ###
		# # ### Appendix A, Step 6(b): update the stick-breaking process 
		# # ###
		
			# # # Appendix A, Step 6(b(i)): create index set I
			# # I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # # Appendix A, Step 6(b(ii)): update eta
			# # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # Appendix A, Step 6(b(ii)): update pi
		    # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # Appendix A, Step 6(b(iii)): update theta
			# # theta <- rgamma(1,priors$r.theta+J-1,priors$q.theta-sum(log(1-eta[-J])))  

			
		# # ###
	    # # ### Appendix A, Step 6(e): update h(t_s)
	    # # ###

		# # mu <- c(mu[I],mu.tmp)
		# # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # exp(log(pie)+get.h(s[x,],S.tilde[mu,3:4],z[x],lc[x],Sigma,Q))))

		# # ###
		# # ###  Appendix A, Step 7: save samples 		   
		# # ###
		
		# # mu.save[,k] <- S.tilde[h,2]
		# # theta.save[k] <- theta    
		# # sigma.mu.save[k] <- sigma.mu
		# # # sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		# # # alpha.save[k,] <- alpha
		# # beta.save[k,] <- beta
		# # gamma.save[k,] <- gamma
		# # v.save[,k] <- v
		# # z.save[,k] <- z
		# # m.save[k] <- m
	# # }
  	
  	# # sigma.mu.save <- sigma.mu.save*s.sd
  	# # t.mcmc.end <- Sys.time()

	# # #####################################################################
	# # ### Write output
	# # #####################################################################
	  
	# # keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# # keep$mu <- keep$mu/sum(m.save)
	# # keep$gamma <- keep$gamma/n.mcmc
	# # cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# # cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# # cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# # cat(paste("\n\nEnd time:",Sys.time()))
	# # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# # cat(paste("\nTime per MCMC iteration:",
		# # round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# # cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	# # list(beta=beta.save,gamma=gamma.save,  # alpha=alpha.save,
		# # mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# # sigma.mu=sigma.mu.save,  # sigma.alpha=sigma.alpha.save,
		# # keep=keep,n.mcmc=n.mcmc)
# # }







# # haulouts.mcmc <- function(s,lc,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,
	# # adapt=FALSE,n.mcmc,n.cores=NULL){
 
 	# # ###
 	# # ### Brian M. Brost (28 DEC 2015)
 	# # ###
 	
 	# # ###
 	# # ### Function arguments: s=telemetry locations; lc=Argos location quality class; 
 	# # ### y=ancillary data source containing binary wet/dry status;
 	# # ### X=design matrix containing covariates influencing wet/dry
 	# # ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# # ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# # ### the location of haul-out sites; S.tilde=matrix summarizing support of haul-out sites 
 	# # ###
 	
 	# # ### See Appendix A of haul-outs manuscript for write-up of this model

	# # ### This algorithm has been modified such that s and y are temporally aligned, and thus
	# # ### z=y. Modifications in the code are commented out and indicated by a "*", i.e., "#*"	

	# # t.start <- Sys.time()
	# # cat(paste("\n\nStart time:",t.start,"\n"))

	# # #####################################################################
	# # ### Libraries and Subroutines
	# # #####################################################################
# # # browser()
	# # truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# # flow <- pnorm(low,mu,sqrt(sig2)) 
		# # fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# # u <- runif(nsamp) 
		# # tmp <- flow+u*(fhigh-flow)
		# # x <- qnorm(tmp,mu,sqrt(sig2))
		# # x
	# # }
	
	# # tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# # (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# # }

	# # get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# # # a <- min(0.01,1/sqrt(k))
		# # a <- min(0.025,1/sqrt(k))
		# # exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# # }

	# # get.Sigma <- function(sigma2,a,rho){  # Get covariance matrix
		# # S <- sigma2*matrix(c(1,sqrt(a)*rho,sqrt(a)*rho,a),2)  # variance-covariance matrix
		# # b <-  S[1,1]*S[2,2]-S[1,2]*S[2,1]  # determinant
		# # P <- (1/b)*matrix(c(S[2,2],-S[2,1],-S[1,2],S[1,1]),2)  # precision matrix
		# # list(P=P,b=b,S=S)
	# # }

	# # dt2 <- function(x,y,z,S,Q,lc=NULL,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # # Density of mixture t-distribution
		# # x <- matrix(x,,2)	
		# # y <- matrix(y,,2)
		# # if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# # if(!is.null(lc)){
			# # S <- S[[lc]]
			# # Q <- Q[[lc]]
		# # }
		# # P <- ifelse(z==1,S,Q)[[1]]  # precision matrix
		# # b <- ifelse(z==1,S$b,Q$b)  # determinant
		# # d <- x-y
		# # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# # if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# # if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# # out
	# # }

	# # # Test dmvt2 function
	# # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # # test2 <- log(0.5)+
		# # # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # plot(test1,test2)
	# # # summary(test2-test1)
	# # # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # # test2 <- 0.5*
		# # # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # plot(test1,test2)
	# # # summary(test2/test1)
	# # # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# # get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # # Accept/reject proposals for mu
		# # mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# # mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# # idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# # idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# # num.z1 <- denom.z1 <- num.z0 <- denom.z0 <- 0
		# # if(length(idx.1)>0){
			# # num.z1 <- sum(sapply(idx.1,function(x)
				# # dt2(s[x,],mu.xy.star,z=1,Sigma,Q,lc[x],log=TRUE)))
			# # denom.z1 <- sum(sapply(idx.1,function(x)
				# # dt2(s[x,],mu.xy,z=1,Sigma,Q,lc[x],log=TRUE)))
		# # }
		# # if(length(idx.0)>0){
			# # num.z0 <- sum(sapply(idx.0,function(x) 
				# # dt2(s[x,],mu.xy.star,z=0,Sigma,Q,lc[x],log=TRUE)))
			# # denom.z0 <- sum(sapply(idx.0,function(x)
				# # dt2(s[x,],mu.xy,z=0,Sigma,Q,lc[x],log=TRUE)))
		# # }
		# # mh.star <- num.z1+num.z0  # +U[mu.star[x],]%*%gamma  # numerator of Metropolis ratio
			# # #-log(sum(exp(U%*%gamma))))  # integral over S.tilde
			# # # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# # mh.0 <-	denom.z1+denom.z0  # +U[mu[x],]%*%gamma   # denominator of Metropolis ratio
			# # #-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# # exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# # }

	# # #####################################################################
	# # ###  Get starting values from previous model
	# # #####################################################################

	# # k <- nrow(start$sigma)
	# # if(!is.null(k)){ 
		# # start$sigma <- start$sigma[k,]
		# # start$a <- start$a[k,]
		# # start$rho <- start$rho[k,]
		# # start$theta <- start$theta[k]+0.1				
		# # start$sigma.mu <- start$sigma.mu[k]
		# # start$gamma <- start$gamma[k,] 
		# # start$h <- start$mu[,k]
		# # start$h <- S.tilde[match(start$h,S.tilde[,2]),1]
		# # start$z <- start$z[,k]
	# # }


	# # #####################################################################
	# # ###  Setup Variables 
	# # #####################################################################
# # # browser() 
	# # cat("\nSetting up variables....")
	# # Ts <- nrow(s)  # number of telemetry locations
	# # #* Ty <- length(y)  # number of wet/dry observations
	# # Ty <- 0
	# # qX <- ncol(X)
	# # #* qW <- ncol(W)
	# # qW <- 0
	# # qU <- ncol(U)
	# # v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# # #* y1 <- which(y==1)+Ts
	# # #* y0 <- which(y==0)+Ts
	# # #* y1.sum <- length(y1)
	# # #* y0.sum <- length(y0)
	# # #* W.cross <- t(W)%*%W  # cross product of W
	# # lc <- as.numeric(lc)  # Argos location quality class
	# # n.lc <- length(unique(lc))  # number of error classes
	# # lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	

	# # #####################################################################
	# # ### Priors
	# # #####################################################################
	
	# # cat("\nGetting priors....")
	
	# # # Observation model
	# # u.sigma <- priors$u.sigma  # upper limit of uniform prior on sigma
	# # mu.sigma <- priors$mu.sigma  # mean of lognormal prior for sigma.mu
	# # sigma.sigma <- priors$sigma.sigma  # variance of longormal prior for sigma.mu
	# # #* r.sigma.alpha <- priors$r.sigma.alpha  # IG prior on sigma.alpha
	# # #* q.sigma.alpha <- priors$q.sigma.alpha  # IG prior on sigma.alpha
	# # r.theta <- priors$r.theta  # IG prior on theta
	# # q.theta <- priors$q.theta  # IG prior on theta
	# # sigma.gamma <- priors$sigma.gamma  # variance of normal prior on gamma
	
	# # # Temporal haul-out process model 
	# # #* mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# # mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# # Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# # Sigma.beta.inv <- solve(Sigma.beta)
	# # A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # # Spatial haul-out process model 
	# # J <- priors$J  # maximum number of clusters per truncation approximation
	# # mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# # # Sigma.gamma <- diag(qU)*sigma.gamma^2
  	# # # Sigma.gamma.inv <- solve(Sigma.gamma)


	# # #####################################################################
	# # ### Standardize parameters
	# # #####################################################################

	# # cat("\nStandardizing variables....")
	
	# # # Center and scale s and S.tilde
	# # s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# # s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# # s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# # S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# # ncol=2,byrow=TRUE))/s.sd
	
	# # # Center and scale tuning parameters
	# # tune$sigma.mu <- tune$sigma.mu/s.sd
	# # tune$mu <- tune$mu/s.sd
	# # tune$sigma <- tune$sigma/s.sd

	# # # Center and scale priors
	# # mu.sigma <- mu.sigma/s.sd  
	# # u.sigma <- u.sigma/s.sd

	# # # Center and scale starting values
	# # sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# # sigma <- start$sigma/s.sd  # observation model standard deviation


	# # #####################################################################
	# # ### Starting values: Appendix A, Steps 1 and 2
	# # #####################################################################

	# # cat("\nGetting starting values....")
# # # browser()
	# # # Observation model
	# # a <- start$a
	# # rho <- start$rho
	# # Sigma <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2,a[x],rho[x]),simplify=FALSE)
	# # Q <- sapply(1:n.lc,function(x) get.Sigma(sigma[x]^2+sigma.mu^2,a[x],rho[x]),simplify=FALSE)

	# # # Temporal haul-out process model
	# # #* alpha <- matrix(start$alpha,qW)  # random effects for temporal haul-out process
	# # #* Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# # #* Sigma.alpha.inv <- solve(Sigma.alpha)
  	# # #* A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# # z <- start$z  # haul-out status; z=1: dry, z=0: wet

	# # # Spatial haul-out process model
	# # gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# # theta <- start$theta  # DP concentration parameter


  	# # #####################################################################
	# # ### Create receptacles for output
	# # #####################################################################
  
  	# # cat("\nCreating receptacles for output....")
	# # sigma.save <- matrix(0,n.mcmc,n.lc)  # longitudianl telemetry measurement error
	# # a.save <- matrix(0,n.mcmc,n.lc)  # adjustment for latitudinal error
	# # rho.save <- matrix(0,n.mcmc,n.lc)  # covariance between long. and lat. errors
	# # beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# # #* alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# # gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# # mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# # theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# # m.save <- numeric(n.mcmc)  # number of clusters
	# # v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# # sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# # #* sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects
	
    
	# # #####################################################################
	# # ### MCMC loop: Appendix A, Steps 3-9
	# # #####################################################################
	
	# # # Track overall MH accpetance rate
	# # keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))

	# # # Track MH accpetance rate for adaptive tuning
	# # keep.tmp <- keep  
	# # T.b <- 50  # frequency of adaptive tuning
	# # m.save.tmp <- 0  # number of clusters

	# # t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# # t.mcmc.start <- Sys.time()  # timing MCMC iterations
	
	# # cat("\nEntering MCMC Loop....\n")
	# # for (k in 1:n.mcmc) {
    	# # if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# # cat(k,"");flush.console()	
    		# # plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# # } 
# # # browser()

		# # if(adapt==TRUE & k%%50==0) {  # Adaptive tuning
			# # # browser()
			# # keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# # keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# # tune$sigma.mu <- get.tune(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# # # tune$gamma <- get.tune(tune$gamma,keep.tmp$gamma,k)
			# # tune$mu <- get.tune(tune$mu,keep.tmp$mu,k)
			# # tune$sigma <- sapply(1:n.lc,function(x) 
				# # get.tune(tune$sigma[x],keep.tmp$sigma[x],k))
			# # tune$a <- sapply(1:n.lc,function(x) get.tune(tune$a[x],keep.tmp$a[x],k))
			# # tune$rho <- sapply(1:n.lc,function(x) get.tune(tune$rho[x],keep.tmp$rho[x],k))
			# # keep.tmp <- lapply(keep.tmp,function(x) x*0)
			# # m.save.tmp <- 0
	   	# # } 	
		
	    
	    # # #--------------------------------------------------------------------------
	  	# # # Update sigma.mu: Appendix A, Step 4 
	    # # #--------------------------------------------------------------------------
# # # browser()
		# # # Lognormal prior
	    # # sigma.mu.star <-  rnorm(1,sigma.mu,tune$sigma.mu)
		# # if(sigma.mu.star>0){
			# # Q.star <- sapply(1:n.lc,function(x) 
				# # get.Sigma(sigma[x]^2+sigma.mu.star^2,a[x],rho[x]),simplify=FALSE)
			# # idx <- which(z==0)
		    # # mh.star.sigma.mu <-	sum(sapply(idx,function(x) 
		    	# # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q.star,lc[x],log=TRUE)))+
		    	# # dnorm(log(sigma.mu.star),log(mu.sigma),sigma.sigma,log=TRUE)
		    # # mh.0.sigma.mu <- sum(sapply(idx,function(x)
		    	# # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=TRUE)))+
		    	# # dnorm(log(sigma.mu),log(mu.sigma),sigma.sigma,log=TRUE)
		    # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# # sigma.mu <- sigma.mu.star
				# # Q <- Q.star
				# # keep$sigma.mu <- keep$sigma.mu+1
				# # keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1
	    	# # } 
		# # }

	
		# # #--------------------------------------------------------------------------
		# # # Update observation model parameters: Appendix A, Step 5
		# # #--------------------------------------------------------------------------
# # # browser()
		# # sigma.star <- rnorm(n.lc,sigma,tune$sigma)  # proposals for sigma
		# # a.star <- rnorm(n.lc,a,tune$a)  # proposals for a
		# # rho.star <- rnorm(n.lc,rho,tune$rho)  # proposals for rho
	
		# # for(i in 1:n.lc){  # loop to iterate over error classes: Appendix A, Step 5(b)

			# # idx <- lc.list[[i]]  # index of locations in error class i
			# # z1 <- idx[which(z[idx]==1)]
			# # z0 <- idx[which(z[idx]==0)]

			# # ###
			# # ### Sample sigma: Appendix A, Step 5(a.i)
			# # ###

			# # if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# # Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])
				# # Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# # mh.star.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# # mh.0.sigma <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# # sigma[i] <- sigma.star[i]
					# # Sigma[[i]] <- Sigma.star
					# # Q[[i]] <- Q.star
					# # keep$sigma[i] <- keep$sigma[i]+1
					# # keep.tmp$sigma[i] <- keep.tmp$sigma[i]+1
				# # }
			# # }

			# # ###
			# # ### Sample a: Appendix A, Step 5(a.ii)
			# # ###
			
			# # if(a.star[i]>0 & a.star[i]<1){
				# # Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])
				# # Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# # mh.star.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# # mh.0.a <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# # if(exp(mh.star.a-mh.0.a)>runif(1)){
					# # a[i] <- a.star[i]
					# # Sigma[[i]] <- Sigma.star
					# # Q[[i]] <- Q.star
					# # keep$a[i] <- keep$a[i]+1
					# # keep.tmp$a[i] <- keep.tmp$a[i]+1
				# # }
			# # }

			# # ###
			# # ### Sample rho: Appendix A, Step 5(a.iii)
			# # ###
			
			# # if(rho.star[i]>0 & rho.star[i]<1){
				# # Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])
				# # Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# # mh.star.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma.star,Q.star))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma.star,Q.star))
				# # mh.0.rho <- sum(dt2(s[z1,],S.tilde[h[z1],3:4],z=1,Sigma,Q,i))+
					# # sum(dt2(s[z0,],S.tilde[h[z0],3:4],z=0,Sigma,Q,i))
				# # if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# # rho[i] <- rho.star[i]
					# # Sigma[[i]] <- Sigma.star
					# # Q[[i]] <- Q.star
					# # keep$rho[i] <- keep$rho[i]+1
					# # keep.tmp$rho[i] <- keep.tmp$rho[i]+1
				# # }
			# # }
		# # }

	
		# # #--------------------------------------------------------------------------
		# # # Update temporal haul-out process model parameters: Appendix A, Step 6
		# # #--------------------------------------------------------------------------
# # # browser()
	 	# # t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# # ###	
		# # ### Update beta: Appendix A, Step 6(a)
		# # ###
		
		# # #* b <- crossprod(X,(v-W%*%alpha))
		# # b <- crossprod(X,v)
	  	# # beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))

		# # ###
		# # ### Update alpha: Appendix A, Step 6(b)
		# # ###
		
		# # #* b <- crossprod(W,(v-X%*%beta))
		# # #* alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# # ###
		# # ### Update v(t_y): Appendix A, Step 6(c)
		# # ###

		# # linpred <- X%*%beta  #* +W%*%alpha  # update linear predictor 
	  	# # #* v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)
	  	# # #* v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)

		# # ###
		# # ### Update v(t_s): Appendix A, Step 6(d)
		# # ###
			
		# # z1 <- which(z==1)
		# # z0 <- which(z==0)
		# # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	  	# # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
	
	  	# # t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of aux. variable update
	  	
	  	# # ###
	  	# # ### Calculate Prob(z(t_s)==1): Appendix A, Step 6(e)
	  	# # ###
	  	
		# # #* p <- pnorm(linpred[1:Ts,])

	    # # ###
	    # # ### Update z: Appendix A, Step 6(f)
	    # # ###
		
		# # #* p1 <- p*sapply(1:Ts,function(x) 
			# # # dt2(s[x,],S.tilde[h[x],3:4],z=1,Sigma,Q,lc[x],log=FALSE))
		# # #* p2 <- (1-p)*sapply(1:Ts,function(x) 
			# # # dt2(s[x,],S.tilde[h[x],3:4],z=0,Sigma,Q,lc[x],log=FALSE))
		# # #* p <- exp(log(p1)-log(p1+p2))	
		# # #* z <- rbinom(Ts,1,p)

		# # ###
		# # ### Update sigma2.alpha: Appendix A, Step 6(g)
		# # ###
		
		# # #* r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/r.sigma.alpha)
		# # #* q.tmp <- qW/2+q.sigma.alpha
		# # #* sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# # #* diag(Sigma.alpha) <- sigma2.alpha
		# # #* Sigma.alpha.inv <- solve(Sigma.alpha)
		# # #* A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		
					
		# # #--------------------------------------------------------------------------
		# # # Update spatial haul-out process model parameters: Appendix A, Step 7
		# # #--------------------------------------------------------------------------
# # # browser()
		# # ###
		# # ### Tabulate cluster membership: Appendix A, Step 7(a) 
		# # ###

		# # n <- table(h)  
		# # m <- length(n)  # number of clusters
		# # mu <- as.numeric(names(n))  # idx of occupied clusters; mu references row in S.tilde

		# # ###
		# # ### Update the stick-breaking process: Appendix A, Step 7(b)
		# # ###
		
			# # # Create index set I: Appendix A, Step 6(b.i)
			# # I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # # Update eta: Appendix A, Step 6(b.ii)
			# # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # Update pi: Appendix A, Step 6(b.ii)
		    # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # Update theta: Appendix A, Step 6(b.iii)
			# # theta <- rgamma(1,r.theta+J-1,q.theta-sum(log(1-eta[-J])))  

		# # ###
		# # ### Update 'unoccupied' mu: Appendix A, Step 7(c)
		# # ###

		# # # p <- exp(U[-mu,]%*%gamma)
		# # # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu
		# # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE)  # idx of unoccupied mu
		
		
		# # ###
		# # ### Update 'occupied' mu: Appendix A, Step 7(d)		   
		# # ###
# # # browser()			
		# # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# # keep$mu <- keep$mu+sum(mh)
		# # keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# # keep.idx <- dup.idx[mh]
		# # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # ###
	    # # ### Update gamma: Appendix A, Step 7(e)
	    # # ###
# # # browser()	    
		# # # Integral over S.tilde, occupied mu only
		# # # gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# # # mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,sigma.gamma,log=TRUE))+
			# # # sum(n*c(U[mu,]%*%gamma.star)-n*log(sum(exp(U%*%gamma.star))))
			# # # # sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # # # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
		# # # mh.0.gamma <- sum(dnorm(gamma,mu.gamma,sigma.gamma,log=TRUE))+
			# # # sum(n*c(U[mu,]%*%gamma)-n*log(sum(exp(U%*%gamma))))
			# # # # sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # # # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# # # if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # # # gamma <- gamma.star
	        # # # # gamma.int <- gamma.int.star
	        # # # keep$gamma <- keep$gamma+1
	        # # # keep.tmp$gamma <- keep.tmp$gamma+1
    	# # # } 
			
		# # ###
	    # # ### Update h(t_s): Appendix A, Step 7(f)
	    # # ###

		# # mu <- c(mu[I],mu.tmp)
		# # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # exp(log(pie)+dt2(s[x,],S.tilde[mu,3:4],z[x],Sigma,Q,lc[x]))))

		# # #--------------------------------------------------------------------------
		# # # Save samples: Appendix A, Step 8
		# # #--------------------------------------------------------------------------
						
		# # sigma.save[k,] <- sigma
		# # a.save[k,] <- a
		# # rho.save[k,] <- rho
		# # mu.save[,k] <- S.tilde[h,2]
		# # theta.save[k] <- theta    
		# # sigma.mu.save[k] <- sigma.mu
		# # #* sigma.alpha.save[k] <- sigma2.alpha
		# # #* alpha.save[k,] <- alpha
		# # beta.save[k,] <- beta
		# # gamma.save[k,] <- gamma
		# # v.save[,k] <- v
		# # z.save[,k] <- z
		# # m.save[k] <- m
		# # m.save.tmp <- m.save.tmp+m
	# # }  # end of MCMC loop

  	# # t.mcmc.end <- Sys.time()
  	
  	# # tune$sigma.mu <- tune$sigma.mu*s.sd
  	# # tune$sigma <- tune$sigma*s.sd
	# # tune$mu <- tune$mu*s.sd
  	# # sigma.save <- sigma.save*s.sd
  	# # sigma.mu.save <- sigma.mu.save*s.sd
	# # #* sigma.alpha.save <- sqrt(sigma.alpha.save)
	
  	# # #####################################################################
	# # ### Write output
	# # #####################################################################
	  
	# # keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# # keep$mu <- keep$mu/sum(m.save)
	# # keep$gamma <- keep$gamma/n.mcmc
	# # keep$sigma <- keep$sigma/n.mcmc
	# # keep$a <- keep$a/n.mcmc
	# # keep$rho <- keep$rho/n.mcmc

	# # cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# # cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# # cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# # cat("\nsigma acceptance rate:",round(keep$sigma,2)) 
	# # cat("\na acceptance rate:",round(keep$a,2))
	# # cat("\nrho acceptance rate:",round(keep$rho,2)) 

	# # cat(paste("\n\nEnd time:",Sys.time()))
	# # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# # cat(paste("\nTime per MCMC iteration:",
		# # round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# # cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))

	# # list(beta=beta.save,gamma=gamma.save,  #* alpha=alpha.save,
		# # mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# # sigma.mu=sigma.mu.save,  #* sigma.alpha=sigma.alpha.save,
		# # sigma=sigma.save,a=a.save,rho=rho.save,
		# # keep=keep,tune=tune,n.mcmc=n.mcmc)
# # }




# # # haulouts.mcmc <- function(s,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# # # ###
 	# # # ### Brian M. Brost (08 DEC 2015)
 	# # # ### Cleaned up version of haulouts.3.mcmc.R
 	# # # ###
 	
 	# # # ###
 	# # # ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# # # ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# # # ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# # # ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# # # ### the location of haul-out sites; S.tilde=support of haul-out sites 
 	# # # ###

	# # # t.start <- Sys.time()
	# # # cat(paste("Start time:",t.start,"\n"))

	# # # #####################################################################
	# # # ### Libraries and Subroutines
	# # # #####################################################################
  
	# # # # library(MCMCpack)  # for Dirichlet distribution functions
	# # # # library(data.table)  # for tabulating and summing
	# # # # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# # # # library(doParallel)  # for parallel processing
	# # # # library(foreach)  # for parallel processing  
	# # # # library(mvtnorm)  # for multivariate normal density
	# # # # library(msm)  # for truncated normal density
  
  	# # # truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# # # flow <- pnorm(low,mu,sqrt(sig2)) 
		# # # fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# # # u <- runif(nsamp) 
		# # # tmp <- flow+u*(fhigh-flow)
		# # # x <- qnorm(tmp,mu,sqrt(sig2))
		# # # x
	# # # }
	
	# # # tkern <- function(d,P,nu=100){  # kernel of t-distribution
		# # # (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# # # }

	# # # adapt <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# # # a <- min(0.01,1/sqrt(k))
		# # # exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	# # # }
	
	# # # dt2 <- function(x,y,z,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
		# # # # browser()		
		# # # x <- matrix(x,,2)	
		# # # y <- matrix(y,,2)
		# # # if(nrow(x)!=nrow(y)) x <- matrix(x,nrow(y),2,byrow=TRUE)
		# # # P <- ifelse(z==1,Sigma,Q)[[1]]  # precision matrix
		# # # b <- ifelse(z==1,Sigma$b,Q$b)  # determinant
		# # # d <- x-y
		# # # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)  # mixture kernel
		# # # if(log) out <- log(out)+log(0.5)+log(b^(-0.5))  # log density
			# # # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# # # if(!log) out <- 0.5*out*b^(-0.5)  # density
			# # # # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
		# # # out
	# # # }

	# # # get.Sigma <- function(sigma2,a,rho){  # get var-cov matrix, determinant, and inverse
		# # # n.lc <- length(sigma2)
		# # # Sigma <- lapply(1:n.lc,function(x) sigma2[x]*  # variance-covariance matrix			
			# # # matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))  
		# # # det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		# # # P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			# # # matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		# # # list(Sigma=Sigma,P=P,det=det)  # return list of lists
	# # # }

	
	
	
	# # # dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # # # Calculate density of mixture t distribution
		# # # x <- matrix(x,,2,byrow=FALSE)
		# # # if(nrow(x)==0) out <- 0
		# # # if(nrow(x)!=0){
			# # # if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# # # lc.idx <- sort(unique(lc))
			# # # n <- length(lc.idx)
			# # # lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# # # P <- Sigma$P[lc.idx]  # precision matrix
			# # # b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))  # determinant
			# # # d <- x-y  
			# # # out <- numeric(nrow(d))
			# # # for(i in 1:n){  # calculate kernel for each Sigma
				# # # idx <- lc.list[[i]]
				# # # d.tmp <- d[idx,]
				# # # out[idx] <- tkern(d.tmp,P[[i]],nu)+tkern(d.tmp,K%*%P[[i]]%*%t(K),nu)
			# # # }	
			# # # if(log) out <- log(0.5)+log(out)+log(b[lc])  # log density
					# # # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			# # # if(!log){ out <- 0.5*out*b[lc]  # density
					# # # # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
			# # # } 
		# # # }
		# # # out
	# # # }

	# # # # Test dmvt2 function
	# # # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # # # test2 <- log(0.5)+
		# # # # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # # plot(test1,test2)
	# # # # summary(test2-test1)
	# # # # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # # # test2 <- 0.5*
		# # # # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # # plot(test1,test2)
	# # # # summary(test2/test1)
	# # # # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# # # get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # # # browser()
		# # # # Accept/reject proposals for mu
		# # # mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# # # mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# # # idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# # # idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# # # mh.star <- sum(  # numerator of Metropolis-Hastings ratio
			# # # dmvt2(s[idx.1,],mu.xy.star,lc[idx.1],Sigma,log=TRUE),
			# # # dmvt2(s[idx.0,],mu.xy.star,lc[idx.0],Q,log=TRUE))+
			# # # U[mu[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# # # mh.0 <-	sum(  # denominator of Metropolis-Hastings ratio
			# # # dmvt2(s[idx.1,],mu.xy,lc[idx.1],Sigma,log=TRUE),
			# # # dmvt2(s[idx.0,],mu.xy,lc[idx.0],Q,log=TRUE))+
			# # # U[mu.star[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# # # exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# # # }

	# # # get.h <- function(x,y,z,lc,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2)){
		# # # # browser()
		# # # # For sampling h
		# # # if(z==1) P <- Sigma$P[[lc]]  # precision matrix when z=1 
		# # # if(z==0) P <- Q$P[[lc]]  # precision matrix when z=0
		# # # d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# # # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# # # log(out)  # log of mixture kernel
	# # # }


	# # # #####################################################################
	# # # ###  Setup Variables 
	# # # #####################################################################
# # # # browser() 
	# # # cat("\nSetting up variables....")

	# # # Ts <- nrow(s)  # number of telemetry locations
	# # # # Ty <- length(y)  # number of wet/dry observations
# # # Ty <- 0
	# # # X <- as.matrix(X)
	# # # U <- as.matrix(U)
	# # # qX <- ncol(X)
	# # # # qW <- ncol(W)
# # # qW <- 0
	# # # # qU <- dim(U)[3]+1	# number of raster layers plus an intercept
# # # qU <- ncol(U)
	# # # v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# # # # y1 <- which(y==1)+Ts
	# # # # y0 <- which(y==0)+Ts
	# # # # y1.sum <- length(y1)
	# # # # y0.sum <- length(y0)
	# # # # W.cross <- t(W)%*%W  # cross product of W
	# # # # idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	# # # # S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# # # # information in S.tilde; note that mu below references row idx in S.tilde
	# # # h <- match(start$h,S.tilde[,2])  # Note: h corresponds to row idx of S.tilde 
		# # # # for computational efficiency, and not idx of mu as in Appendix A
	# # # # U <- cbind(1,values(U)[idx])  # convert raster to design matrix

	# # # lc <- as.numeric(priors$lc)  # Argos location quality class
	# # # n.lc <- length(unique(lc))  # number of error classes
	# # # lc.list <- sapply(sort(unique(lc)),function(x) which(lc==x),simplify=FALSE)
	
	# # # #####################################################################
	# # # ### Standardize parameters
	# # # #####################################################################

	# # # cat("\nStandardizing variables....")
	
	# # # # Center and scale s and S.tilde
	# # # s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# # # s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# # # s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# # # S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# # # ncol=2,byrow=TRUE))/s.sd
	
	# # # # Center and scale tuning parameters
	# # # tune$sigma.mu <- tune$sigma.mu/s.sd
	# # # tune$mu <- tune$mu/s.sd
	# # # tune$sigma <- tune$sigma/s.sd

	# # # # Center and scale priors
	# # # priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	# # # priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu
	# # # priors$mu.sigma <- priors$mu.sigma/s.sd  # variance on lognormal prior for sigma.mu

	# # # # Center and scale starting values
	# # # sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter
	# # # sigma <- start$sigma/s.sd  # observation model standard deviation


	# # # #####################################################################
	# # # ### Priors
	# # # #####################################################################
	
	# # # cat("\nGetting priors....")
	
	# # # # Observation model
	# # # u.sigma <- priors$u.sigma/s.sd
	
	# # # # Temporal haul-out process model 
	# # # # mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# # # mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# # # Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# # # Sigma.beta.inv <- solve(Sigma.beta)
	# # # A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # # # Spatial haul-out process model 
	# # # J <- priors$J  # maximum number of clusters per truncation approximation
	# # # mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# # # Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	# # # Sigma.gamma.inv <- solve(Sigma.gamma)


	# # # #####################################################################
	# # # ### Appendix A, Step 1: starting values 
	# # # #####################################################################

	# # # cat("\nGetting starting values....")

	# # # # Observation model
	# # # a <- start$a
	# # # rho <- start$rho
	# # # Sigma <- get.Sigma(sigma^2,a,rho)
	# # # Q <- get.Sigma(sigma^2+sigma.mu^2,a,rho)

	# # # # Temporal haul-out process model
	# # # beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	# # # # alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	# # # # Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# # # # Sigma.alpha.inv <- solve(Sigma.alpha)
  	# # # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# # # z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # # # linpred <- X%*%beta+W%*%alpha
# # # linpred <- X%*%beta

	# # # # Spatial haul-out process model
	# # # gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # # # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# # # theta <- start$theta  # DP concentration parameter


  	# # # #####################################################################
	# # # ### Create receptacles for output
	# # # #####################################################################
  
  	# # # cat("\nCreating receptacles for output....")
	# # # beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# # # # alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# # # gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# # # mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# # # theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# # # m.save <- numeric(n.mcmc)  # number of clusters
	# # # v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# # # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# # # sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# # # # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

	# # # sigma.save <- matrix(0,n.mcmc,n.lc)
	# # # a.save <- matrix(0,n.mcmc,n.lc)
	# # # rho.save <- matrix(0,n.mcmc,n.lc)

    
	# # # #####################################################################
	# # # ### Appendix A, Steps 3-8: MCMC loop 
	# # # #####################################################################

	# # # keep <- list(mu=0,sigma.mu=0,gamma=0)  # number of proposals accepted for Metropolis updates

# # # keep <- list(mu=0,sigma.mu=0,gamma=0,sigma=rep(0,n.lc),a=rep(0,n.lc),rho=rep(0,n.lc))
	# # # keep.tmp <- list(mu=0,sigma.mu=0,gamma=0)  # for adaptive tuning
	# # # m.save.tmp <- 0
	# # # t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# # # t.mcmc.start <- Sys.time()  # timing MCMC iterations
	# # # T.b <- 50  # frequency of adaptive tuning
	
	# # # cat("\nEntering MCMC Loop....\n")
	# # # for (k in 1:n.mcmc) {
    	# # # if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# # # cat(k,"");flush.console()	
    		# # # plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# # # } 

    	# # # if(k%%50==0) {  # Adaptive tuning
			# # # # browser()
			# # # keep.tmp$mu <- keep.tmp$mu/m.save.tmp
			# # # keep.tmp[-1] <- lapply(keep.tmp[-1],function(x) x/T.b)
			# # # tune$sigma.mu <- adapt(tune$sigma.mu,keep.tmp$sigma.mu,k)
			# # # tune$gamma <- adapt(tune$gamma,keep.tmp$gamma,k)
			# # # tune$mu <- adapt(tune$mu,keep.tmp$mu,k)
			# # # keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	# # # } 
	
		# # # #--------------------------------------------------------------------------
		# # # # Appendix A, Step 4: update observation model parameters (Sigma)
		# # # #--------------------------------------------------------------------------

# # # # browser()
# # # # priors$sigma <- priors$sigma/s.scale

		# # # sigma.star <- rnorm(n.lc,sigma,tune$sigma) #Proposals for sigma
		# # # a.star <- rnorm(n.lc,a,tune$a) #Proposals for a
		# # # rho.star <- rnorm(n.lc,rho,tune$rho) #Proposals for rho

	# # # dmvt2.vec <- function(d,P,b,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # # # Calculate density of mixture t distribution
		# # # b <- b^(-0.5)  # determinant
		# # # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# # # out <- log(0.5)+log(out)+log(b)  # log density
					# # # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
		# # # out
	# # # }
	
		# # # for(i in 1:n.lc){ #Loop to iterate over error classes: Appendix A, step 2(f)
	# # # # i <- 2
			# # # idx <- lc.list[[i]] #Index of locations in error class i
			# # # d.tmp <- s[idx,]-S.tilde[h[idx],3:4]
			# # # K <- matrix(c(-1,0,0,1),2)
			# # # z1 <- which(z[idx]==1)
			# # # z0 <- which(z[idx]==0)
# # # # browser()	
			# # # ### Sample sigma: Appendix A, step 2(b)

			# # # if(sigma.star[i]>0 & sigma.star[i]<u.sigma){
				# # # Q.star <- get.Sigma(sigma.star[i]^2+sigma.mu^2,a[i],rho[i])
				# # # Sigma.star <- get.Sigma(sigma.star[i]^2,a[i],rho[i])

				# # # mh.star.sigma <- 
				# # # sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# # # mh.0.sigma <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))

		# # # # mh.star.sigma <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		# # # # mh.0.sigma <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))

				# # # if(exp(mh.star.sigma-mh.0.sigma)>runif(1)){
					# # # sigma[i] <- sigma.star[i]
					# # # keep$sigma[i] <- keep$sigma[i]+1
					# # # Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# # # Sigma$P[[i]] <- Sigma.star$P[[1]]
					# # # Sigma$det[[i]] <- Sigma.star$det[[1]]
					# # # Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# # # Q$P[[i]] <- Q.star$P[[1]]
					# # # Q$det[[i]] <- Q.star$det[[1]]
				# # # }
			# # # }

			# # # ### Sample a: Appendix A, step 2(c)

			# # # if(a.star[i]>0 & a.star[i]<1){
				# # # Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a.star[i],rho[i])
				# # # Sigma.star <- get.Sigma(sigma[i]^2,a.star[i],rho[i])

				# # # mh.star.a <- 
				# # # sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# # # mh.0.a <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))
				# # # # mh.star.a <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
			  	# # # # mh.0.a <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
				# # # if(exp(mh.star.a-mh.0.a)>runif(1)){
					# # # a[i] <- a.star[i]
					# # # keep$a[i] <- keep$a[i]+1
					# # # Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# # # Sigma$P[[i]] <- Sigma.star$P[[1]]
					# # # Sigma$det[[i]] <- Sigma.star$det[[1]]
					# # # Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# # # Q$P[[i]] <- Q.star$P[[1]]
					# # # Q$det[[i]] <- Q.star$det[[1]]
				# # # }
			# # # }

			# # # ### Sample rho: Appendix A, step 2(d)

			# # # if(rho.star[i]>0 & rho.star[i]<1){
				# # # Q.star <- get.Sigma(sigma[i]^2+sigma.mu^2,a[i],rho.star[i])
				# # # Sigma.star <- get.Sigma(sigma[i]^2,a[i],rho.star[i])

				# # # mh.star.rho <- 
				# # # sum(dmvt2.vec(d.tmp[z1,],Sigma.star$P[[1]],b=Sigma.star$det[[1]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q.star$P[[1]],b=Q.star$det[[1]]))
				# # # mh.0.rho <- sum(dmvt2.vec(d.tmp[z1,],Sigma$P[[i]],b=Sigma$det[[i]]))+
				# # # sum(dmvt2.vec(d.tmp[z0,],Q$P[[i]],b=Q$det[[i]]))

				# # # # mh.star.rho <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
			  	# # # # mh.0.rho <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
				# # # if(exp(mh.star.rho-mh.0.rho)>runif(1)){
					# # # rho[i] <- rho.star[i]
					# # # keep$rho[i] <- keep$rho[i]+1
					# # # Sigma$Sigma[[i]] <- Sigma.star$Sigma[[1]]
					# # # Sigma$P[[i]] <- Sigma.star$P[[1]]
					# # # Sigma$det[[i]] <- Sigma.star$det[[1]]
					# # # Q$Sigma[[i]] <- Q.star$Sigma[[1]]
					# # # Q$P[[i]] <- Q.star$P[[1]]
					# # # Q$det[[i]] <- Q.star$det[[1]]
				# # # }
			# # # }
		# # # }


	
		# # # #--------------------------------------------------------------------------
		# # # # Appendix A, Step 4: update temporal haul-out process model parameters 
		# # # #--------------------------------------------------------------------------
# # # # browser()
	 	# # # t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# # # ###
		# # # ### Appendix A, Step 4(a): update v(t_y)
		# # # ###
		
	  	# # # # v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# # # # v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# # # ###
		# # # ### Appendix A, Step 4(b): update v(t_s)
		# # # ###
			
		# # # z1 <- which(z==1)
		# # # z0 <- which(z==0)
	  	# # # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# # # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	
		# # # ###
		# # # ### Appendix A, Step 4(c): update alpha
		# # # ###
		
		# # # # b <- crossprod(W,(v-X%*%beta))
		# # # # alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# # # ###	
		# # # ### Appendix A, Step 4(d): update beta
		# # # ###
		
		# # # # b <- crossprod(X,(v-W%*%alpha))
# # # # browser()
# # # b <- crossprod(X,v)
	  	# # # beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	# # # t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of auxilliary variable update
	  	
	  	# # # ###
	  	# # # ### Appendix A, Step 4(e): calculate Prob(z(t_s)==1)
	  	# # # ###
	  	
	  	# # # # linpred <- X%*%beta+W%*%alpha  # update linear predictor 
		# # # # p <- pnorm(linpred[1:Ts,])

	    # # # ###
	    # # # ### Appendix A, Step 4(f): update z
	    # # # ###
	    
		# # # # p1 <- p*dmvt2(s,S.tilde[h,3:4],lc,Sigma,log=FALSE)
		# # # # p2 <- (1-p)*dmvt2(s,S.tilde[h,3:4],lc,Q,log=FALSE)
		# # # # p <- exp(log(p1)-log(p1+p2))	
		# # # # z <- rbinom(Ts,1,p)

		# # # ###
		# # # ### Appendix A, Step 4(g): update sigma2.alpha
		# # # ###
		
		# # # # r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		# # # # q.tmp <- qW/2+priors$q.sigma.alpha
		# # # # sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# # # # diag(Sigma.alpha) <- sigma2.alpha
		# # # # Sigma.alpha.inv <- solve(Sigma.alpha)
		# # # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		

	    # # # #--------------------------------------------------------------------------
	  	# # # # Appendix A, Step 5: update sigma.mu
	    # # # #--------------------------------------------------------------------------
# # # # browser()

		# # # # test <- rnorm(1000,log(mu.sigma),tau)
		# # # # hist(s.sd*exp(test))

		# # # # Lognormal prior
	    # # # sigma.mu.star <-  exp(rnorm(1,log(sigma.mu),tune$sigma.mu))
		# # # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,a,rho)
		# # # idx <- which(z==0)
	    # # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))					+dnorm(log(sigma.mu.star),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))						    +dnorm(log(sigma.mu),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	# # # sigma.mu <- sigma.mu.star
			# # # Q <- Q.star
			# # # keep$sigma.mu <- keep$sigma.mu+1
			# # # keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1

    	# # # } 

		# # # # Uniform prior	    
	    # # # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)				   
	    # # # # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# # # # # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,n.lc,Mix)
			# # # # Q.star <- get.Sigma(sigma^2+sigma.mu.star^2,a,rho)
			# # # # idx <- which(z==0)
		    # # # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		    # # # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
		    # # # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# # # # sigma.mu <- sigma.mu.star
				# # # # Q <- Q.star
				# # # # keep$sigma.mu <- keep$sigma.mu+1
	    	# # # # } 
	    # # # # }
		
				
		# # # #--------------------------------------------------------------------------
		# # # # Update spatial haul-out process model parameters: Appendix A, Step 6
		# # # #--------------------------------------------------------------------------
# # # # browser()
		# # # ###
		# # # ### Appendix A, Step 6(a): tabulate cluster membership 
		# # # ###

		# # # n <- table(h)  
		# # # m <- length(n)  # number of clusters
		# # # mu <- as.numeric(names(n))  # idx of occupied clusters

		# # # ###
		# # # ### Appendix A, Step 6(d): update 'unoccupied' mu
		# # # ###
# # # # browser()
		# # # p <- exp(U[-mu,]%*%gamma)
		# # # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		# # # ###
		# # # ### Appendix A, Step 6(c): update 'occupied' mu		   
		# # # ###
# # # # browser()			
		# # # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# # # keep$mu <- keep$mu+sum(mh)
		# # # keep.tmp$mu <- keep.tmp$mu+sum(mh)
		# # # keep.idx <- dup.idx[mh]
		# # # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # # ###
	    # # # ### Appendix A, Step 6(f): update gamma
	    # # # ###

		# # # # Integral over S.tilde, occupied mu only
		# # # gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# # # mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # # sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # # # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
			# # # mh.0.gamma <- sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # # sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # # # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# # # if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # # # gamma <- gamma.star
	        # # # # gamma.int <- gamma.int.star
	        # # # keep$gamma <- keep$gamma+1
	        # # # keep.tmp$gamma <- keep.tmp$gamma+1
    	# # # } 

		# # # ###
		# # # ### Appendix A, Step 6(b): update the stick-breaking process 
		# # # ###
		
			# # # # Appendix A, Step 6(b(i)): create index set I
			# # # I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # # # Appendix A, Step 6(b(ii)): update eta
			# # # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# # # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # # Appendix A, Step 6(b(ii)): update pi
		    # # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # # Appendix A, Step 6(b(iii)): update theta
			# # # theta <- rgamma(1,priors$r.theta+J-1,priors$q.theta-sum(log(1-eta[-J])))  

			
		# # # ###
	    # # # ### Appendix A, Step 6(e): update h(t_s)
	    # # # ###

		# # # mu <- c(mu[I],mu.tmp)
		# # # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # # exp(log(pie)+get.h(s[x,],S.tilde[mu,3:4],z[x],lc[x],Sigma,Q))))

		# # # ###
		# # # ###  Appendix A, Step 7: save samples 		   
		# # # ###
		
		# # # mu.save[,k] <- S.tilde[h,2]
		# # # theta.save[k] <- theta    
		# # # sigma.mu.save[k] <- sigma.mu
		# # # # sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		# # # # alpha.save[k,] <- alpha
		# # # beta.save[k,] <- beta
		# # # gamma.save[k,] <- gamma
		# # # v.save[,k] <- v
		# # # z.save[,k] <- z
		# # # m.save[k] <- m
		# # # m.save.tmp <- m.save.tmp+m
# # # sigma.save[k,] <- sigma*s.sd
# # # a.save[k,] <- a
# # # rho.save[k,] <- rho
	# # # }
  	
  	# # # tune$sigma.mu <- tune$sigma.mu*s.sd
	# # # tune$mu <- tune$mu*s.sd
  	# # # sigma.mu.save <- sigma.mu.save*s.sd
  	# # # t.mcmc.end <- Sys.time()

	# # # #####################################################################
	# # # ### Write output
	# # # #####################################################################
	  
	# # # keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# # # keep$mu <- keep$mu/sum(m.save)
	# # # keep$gamma <- keep$gamma/n.mcmc
	# # # keep$sigma <- keep$sigma/n.mcmc
	# # # keep$a <- keep$a/n.mcmc
	# # # keep$rho <- keep$rho/n.mcmc

	# # # cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# # # cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# # # cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 

	# # # cat(paste("\nsigma acceptance rate:",round(keep$sigma,2))) 
	# # # cat(paste("\na acceptance rate:",round(keep$a,2))) 
	# # # cat(paste("\nrho acceptance rate:",round(keep$rho,2))) 

	# # # cat(paste("\n\nEnd time:",Sys.time()))
	# # # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# # # cat(paste("\nTime per MCMC iteration:",
		# # # round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# # # cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	# # # list(beta=beta.save,gamma=gamma.save,  # alpha=alpha.save,
		# # # mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# # # sigma.mu=sigma.mu.save,  # sigma.alpha=sigma.alpha.save,
		# # # keep=keep,tune=tune,n.mcmc=n.mcmc,
		# # # sigma=sigma.save,a=a.save,rho=rho.save)
# # # }



# # # haulouts.mcmc <- function(s,y=NULL,X,W=NULL,U,S.tilde,priors,tune,start,n.mcmc,n.cores=NULL){
 
 	# # # ###
 	# # # ### Brian M. Brost (08 DEC 2015)
 	# # # ### See haulouts.sim.R to simulate data according to this model specification,
 	# # # ### and haulouts.script.R for implementations with harbor seal data
 	# # # ###
 	
 	# # # ###
 	# # # ### Function arguments: s=telemetry locations, y=ancillary data source containing
 	# # # ### binary wet/dry status; X=design matrix containing covariates influencing wet/dry
 	# # # ### status of telemetry locations s and wet/dry status of ancillary data y;  
 	# # # ### W=basis expansion for s and y; U=design matrix containing covariates influencing
 	# # # ### the location of haul-out sites; S.tilde=support of haul-out sites 
 	# # # ###

	# # # t.start <- Sys.time()
	# # # cat(paste("Start time:",t.start,"\n"))

	# # # #####################################################################
	# # # ### Libraries and Subroutines
	# # # #####################################################################
  
	# # # # library(MCMCpack)  # for Dirichlet distribution functions
	# # # # library(data.table)  # for tabulating and summing
	# # # # library(dplyr)  # dense_rank() for ranking clusters smallest to largest
	# # # # library(doParallel)  # for parallel processing
	# # # # library(foreach)  # for parallel processing  
	# # # # library(mvtnorm)  # for multivariate normal density
	# # # # library(msm)  # for truncated normal density
  
  	# # # truncnormsamp <- function(mu,sig2,low,high,nsamp){  # truncated normal sampler
		# # # flow <- pnorm(low,mu,sqrt(sig2)) 
		# # # fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# # # u <- runif(nsamp) 
		# # # tmp <- flow+u*(fhigh-flow)
		# # # x <- qnorm(tmp,mu,sqrt(sig2))
		# # # x
	# # # }
	
	# # # tkern <- function(d,P,nu){  # kernel of t-distribution
		# # # (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2)
	# # # }

	# # # get.Sigma <- function(sigma2,a,rho){  # get var-cov matrix, determinant, and inverse
		# # # n.lc <- length(sigma2)
		# # # Sigma <- lapply(1:n.lc,function(x) sigma2[x]*  # variance-covariance matrix			
			# # # matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))  
		# # # det <- lapply(Sigma,function(x) x[1,1]*x[2,2]-x[1,2]*x[2,1])  # determinant
		# # # P <- lapply(1:n.lc,function(x) (1/det[[x]])*  # precision matrix
			# # # matrix(c(Sigma[[x]][2,2],-Sigma[[x]][2,1],-Sigma[[x]][1,2],Sigma[[x]][1,1]),2))
		# # # list(Sigma=Sigma,P=P,det=det)  # return list of lists
	# # # }

	# # # dmvt2 <- function(x,y,lc,Sigma,nu=100,K=matrix(c(-1,0,0,1),2),log=TRUE){
	 	# # # # Calculate density of mixture t distribution
		# # # x <- matrix(x,,2,byrow=FALSE)
		# # # if(nrow(x)==0) out <- 0
		# # # if(nrow(x)!=0){
			# # # if(!is.matrix(y)) y <- matrix(y,nrow(x),2,byrow=TRUE)
			# # # lc.idx <- sort(unique(lc))
			# # # n <- length(lc.idx)
			# # # lc.list <- sapply(lc.idx,function(x) which(lc==x),simplify=FALSE)
			# # # P <- Sigma$P[lc.idx]  # precision matrix
			# # # b <- unlist(lapply(Sigma$det,function(x) x^(-0.5)))  # determinant
			# # # d <- x-y  
			# # # out <- numeric(nrow(d))
			# # # for(i in 1:n){  # calculate kernel for each Sigma
				# # # idx <- lc.list[[i]]
				# # # d.tmp <- d[idx,]
				# # # out[idx] <- tkern(d.tmp,P[[i]],nu)+tkern(d.tmp,K%*%P[[i]]%*%t(K),nu)
				# # # # out[idx] <- (1+1/nu*(rowSums((d.tmp%*%P[[i]])*d.tmp)))^-((nu+2)/2) +
					# # # # (1+1/nu*(rowSums((d.tmp%*%(K%*%P[[i]]%*%t(K)))*d.tmp)))^-((nu+2)/2)
			# # # }	
			# # # if(log) out <- log(0.5)+log(out)+log(b[lc])  # log density
					# # # # +lgamma((nu+2)/2)-(lgamma(nu/2)+log(nu)+log(pi))  # constant
			# # # if(!log){ out <- 0.5*out*b[lc]  # density
					# # # # *gamma((nu+2)/2)/(gamma(nu/2)*nu*pi)  # constant
			# # # } 
		# # # }
		# # # out
	# # # }

	# # # # Test dmvt2 function
	# # # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=TRUE)
	# # # # test2 <- log(0.5)+
		# # # # log(sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # # plot(test1,test2)
	# # # # summary(test2-test1)
	# # # # lgamma((100+2)/2)-(lgamma(100/2)+log(100)+log(pi))
	
	# # # # test1 <- dmvt2(s,S.tilde[ht,3:4],lc,Sigma,nu=100,K=K,log=FALSE)
	# # # # test2 <- 0.5*
		# # # # (sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],sigma=Sigma$Sigma[[lc[x]]],
		# # # # df=100,log=FALSE))+sapply(1:Ts,function(x) dmvt(s[x,],S.tilde[ht[x],3:4],
		# # # # sigma=K%*%Sigma$Sigma[[lc[x]]]%*%t(K),df=100,log=FALSE)))
	# # # # plot(test1,test2)
	# # # # summary(test2/test1)
	# # # # gamma((100+2)/2)/(gamma(100/2)*100*pi)
  
	# # # get.mh.mu <- function(x,mu,mu.star,mu.tmp,S,h,s,lc,z,Sigma,Q,U,gamma){
		# # # # browser()
		# # # # Accept/reject proposals for mu
		# # # mu.xy <- S[mu[x],3:4]  # location of current clusters mu
		# # # mu.xy.star <- S[mu.star[x],3:4]  # location of proposal clusters mu.star
		# # # idx.0 <- which(h==mu[x]&z==0)  # obs. associated with mu and z=0
		# # # idx.1 <- which(h==mu[x]&z==1)  # obs. associated with mu and z=1
		# # # mh.star <- sum(  # numerator of Metropolis-Hastings ratio
			# # # dmvt2(s[idx.1,],mu.xy.star,lc[idx.1],Sigma,log=TRUE),
			# # # dmvt2(s[idx.0,],mu.xy.star,lc[idx.0],Q,log=TRUE))+
			# # # U[mu[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # # log(sum(exp(U[c(mu.star[x],mu[-x],mu.tmp),]%*%gamma)))  # integral over active mu
		# # # mh.0 <-	sum(  # denominator of Metropolis-Hastings ratio
			# # # dmvt2(s[idx.1,],mu.xy,lc[idx.1],Sigma,log=TRUE),
			# # # dmvt2(s[idx.0,],mu.xy,lc[idx.0],Q,log=TRUE))+
			# # # U[mu.star[x],]%*%gamma#-log(sum(exp(U%*%gamma)))  # integral over S.tilde
			# # # # log(sum(exp(U[c(mu,mu.tmp),]%*%gamma)))  # integral over active mu
		# # # exp(mh.star-mh.0)>runif(1)  # Accept or reject
	# # # }

	# # # get.h <- function(x,y,z,lc,Sigma,Q,nu=100,K=matrix(c(-1,0,0,1),2)){
		# # # # browser()
		# # # # For sampling h
		# # # if(z==1) P <- Sigma$P[[lc]]  # precision matrix when z=1 
		# # # if(z==0) P <- Q$P[[lc]]  # precision matrix when z=0
		# # # d <- matrix(x,nrow(y),2,byrow=TRUE)-y
		# # # out <- tkern(d,P,nu)+tkern(d,K%*%P%*%t(K),nu)
		# # # # out <- (1+1/nu*(rowSums((d%*%P)*d)))^-((nu+2)/2) +
			# # # # (1+1/nu*(rowSums((d%*%(K%*%P%*%t(K)))*d)))^-((nu+2)/2)
		# # # log(out)  # log of mixture kernel
	# # # }


	# # # #####################################################################
	# # # ###  Setup Variables 
	# # # #####################################################################
# # # # browser() 
	# # # cat("\nSetting up variables....")

	# # # Ts <- nrow(s)  # number of telemetry locations
	# # # # Ty <- length(y)  # number of wet/dry observations
# # # Ty <- 0
	# # # X <- as.matrix(X)
	# # # U <- as.matrix(U)
	# # # qX <- ncol(X)
	# # # # qW <- ncol(W)
# # # qW <- 0
	# # # # qU <- dim(U)[3]+1	# number of raster layers plus an intercept
# # # qU <- ncol(U)
	# # # v <- numeric(Ty+Ts)  # auxilliary variable for continuous haul-out process
	# # # # y1 <- which(y==1)+Ts
	# # # # y0 <- which(y==0)+Ts
	# # # # y1.sum <- length(y1)
	# # # # y0.sum <- length(y0)
	# # # # W.cross <- t(W)%*%W  # cross product of W
	# # # # idx <- which(values(S.tilde)==1)  # cells that define S.tilde
	# # # # S.tilde <- cbind(1:length(idx),idx,xyFromCell(S.tilde,idx))  # matrix summarizing
		# # # # information in S.tilde; note that mu below references row idx in S.tilde
	# # # h <- match(start$h,S.tilde[,2])  # Note: h corresponds to row idx of S.tilde 
		# # # # for computational efficiency, and not idx of mu as in Appendix A
	# # # # U <- cbind(1,values(U)[idx])  # convert raster to design matrix
	
	
	# # # #####################################################################
	# # # ### Standardize parameters
	# # # #####################################################################

	# # # cat("\nStandardizing variables....")
	
	# # # # Center and scale s and S.tilde
	# # # s.sd <- max(apply(s,2,function(x) max(x)-min(x)))/6
	# # # s.mean <- apply(s,2,function(x) max(x)+min(x))/2
	# # # s <- (s-matrix(s.mean,nrow=Ts,ncol=2,byrow=TRUE))/s.sd
	# # # S.tilde[,3:4] <- (S.tilde[,3:4]-matrix(s.mean,nrow=nrow(S.tilde),
		# # # ncol=2,byrow=TRUE))/s.sd
	
	# # # # Center and scale tuning parameters
	# # # tune$sigma.mu <- tune$sigma.mu/s.sd
	# # # tune$mu <- tune$mu/s.sd

	# # # # Center and scale priors
	# # # priors$sigma <- priors$sigma/s.sd  # observation model standard deviation
	# # # priors$sigma.mu.l <- priors$sigma.mu.l/s.sd  # lower bound of uniform prior on sigma.mu
	# # # priors$sigma.mu.u <- priors$sigma.mu.u/s.sd  # upper bound of uniform prior on sigma.mu

	# # # # Center and scale starting values
	# # # sigma.mu <- start$sigma.mu/s.sd  # homerange dispersion parameter


	# # # #####################################################################
	# # # ### Priors
	# # # #####################################################################
	
	# # # cat("\nGetting priors....")
	
	# # # # Observation model
	# # # lc <- as.numeric(priors$lc)  # Argos location quality class
	# # # sigma2 <- priors$sigma^2  # standardized variance component
	# # # Sigma <- get.Sigma(sigma2,priors$a,priors$rho)
	# # # Q <- get.Sigma(sigma2+sigma.mu^2,priors$a,priors$rho)
# # # # browser()
	# # # # Temporal haul-out process model 
	# # # mu.alpha <- matrix(0,qW,1)  # random effects for temporal haul-out process
	# # # mu.beta <- matrix(0,qX,1)  # temporal haul-out process coefficients; fixed effects
	# # # Sigma.beta <- diag(qX)*priors$sigma.beta^2
	# # # Sigma.beta.inv <- solve(Sigma.beta)
	# # # A.inv.beta <- solve(t(X)%*%X+Sigma.beta.inv)

	# # # # Spatial haul-out process model 
	# # # J <- priors$J  # maximum number of clusters per truncation approximation
	# # # mu.gamma <- matrix(0,qU,1)  # haul-out location RSF coefficients
	# # # Sigma.gamma <- diag(qU)*priors$sigma.gamma^2
  	# # # Sigma.gamma.inv <- solve(Sigma.gamma)


	# # # #####################################################################
	# # # ### Appendix A, Step 1: starting values 
	# # # #####################################################################

	# # # cat("\nGetting starting values....")

	# # # # Temporal haul-out process model
	# # # beta <- matrix(start$beta,qX)  # temporal haul-out process coefficients; fixed effects
	# # # # alpha <- matrix(0,qW)  # random effects for temporal haul-out process
	# # # # Sigma.alpha <- diag(qW)*start$sigma.alpha^2
  	# # # # Sigma.alpha.inv <- solve(Sigma.alpha)
  	# # # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
	# # # z <- start$z  # haul-out status; z=1: dry, z=0: wet
	# # # # linpred <- X%*%beta+W%*%alpha
# # # linpred <- X%*%beta
	# # # # Spatial haul-out process model
	# # # gamma <- matrix(start$gamma,qU)  # haul-out location RSF coefficients
	# # # # gamma.int <- log(sum(exp(U%*%gamma)))  # integral for denominator of RSF
	# # # theta <- start$theta  # DP concentration parameter


  	# # # #####################################################################
	# # # ### Create receptacles for output
	# # # #####################################################################
  
  	# # # cat("\nCreating receptacles for output....")
	# # # beta.save <- matrix(0,n.mcmc,qX)  # temporal haul-out process coefficients; fixed effects
	# # # alpha.save <- matrix(0,n.mcmc,qW)  # random effects of temporal haul-out process
	# # # gamma.save <- matrix(0,n.mcmc,qU)  # haul-out location RSF coefficients
	# # # mu.save <- matrix(0,Ts,n.mcmc)  # DP cluster assignment indicator
	# # # theta.save <- numeric(n.mcmc)  # DP concentration parameter
	# # # m.save <- numeric(n.mcmc)  # number of clusters
	# # # v.save <- matrix(0,Ts+Ty,n.mcmc)  # auxiliary variable for temporal haul-out process
	# # # z.save <- matrix(0,Ts,n.mcmc)  # haul-out indicator variable
	# # # sigma.mu.save <- numeric(n.mcmc)  # haul-out dispersion parameter
	# # # sigma.alpha.save <- numeric(n.mcmc)  # standard deviation of random effects

    
	# # # #####################################################################
	# # # ### Appendix A, Steps 3-8: MCMC loop 
	# # # #####################################################################

	# # # keep <- list(mu=0,sigma.mu=0,gamma=0)  # number of proposals accepted for Metropolis updates
	# # # t.v.update <- 0  # timing updates of auxiliary variable for temporal haul-out process
	# # # t.mcmc.start <- Sys.time()  # timing MCMC iterations

  	# # # cat("\nEntering MCMC Loop....\n")
	# # # for (k in 1:n.mcmc) {
    	# # # if(k%%1000==0) {  # Monitor the appropriateness of J, the truncation approximation
	    	# # # cat(k,"");flush.console()	
    		# # # plot(pie,type="b",ylab=expression(pi),xlab="Cluster index",las=1,pch=19,cex=0.5)
    	# # # } 
	
	
		# # # #--------------------------------------------------------------------------
		# # # # Appendix A, Step 4: update temporal haul-out process model parameters 
		# # # #--------------------------------------------------------------------------
# # # # browser()
	 	# # # t.v.start <- Sys.time()  # start time of auxiliary variable update
	
		# # # ###
		# # # ### Appendix A, Step 4(a): update v(t_y)
		# # # ###
		
	  	# # # # v[y1] <- truncnormsamp(linpred[y1],1,0,Inf,y1.sum)
	  	# # # # v[y0] <- truncnormsamp(linpred[y0],1,-Inf,0,y0.sum)

		# # # ###
		# # # ### Appendix A, Step 4(b): update v(t_s)
		# # # ###
			
		# # # z1 <- which(z==1)
		# # # z0 <- which(z==0)
	  	# # # v[z1] <- truncnormsamp(linpred[z1],1,0,Inf,length(z1))
		# # # v[z0] <- truncnormsamp(linpred[z0],1,-Inf,0,length(z0))
	
		# # # ###
		# # # ### Appendix A, Step 4(c): update alpha
		# # # ###
		
		# # # # b <- crossprod(W,(v-X%*%beta))
		# # # # alpha <- A.inv.alpha%*%b+crossprod(chol(A.inv.alpha),matrix(rnorm(qW),qW,1))

		# # # ###	
		# # # ### Appendix A, Step 4(d): update beta
		# # # ###
		
		# # # # b <- crossprod(X,(v-W%*%alpha))
# # # b <- crossprod(X,v)
	  	# # # beta <- A.inv.beta%*%b+crossprod(chol(A.inv.beta),matrix(rnorm(qX),qX,1))
	  	
	  	# # # t.v.update <- t.v.update+difftime(Sys.time(),t.v.start,"secs")  # end time of auxilliary variable update
	  	
	  	# # # ###
	  	# # # ### Appendix A, Step 4(e): calculate Prob(z(t_s)==1)
	  	# # # ###
	  	
	  	# # # # linpred <- X%*%beta+W%*%alpha  # update linear predictor 
		# # # # p <- pnorm(linpred[1:Ts,])

	    # # # ###
	    # # # ### Appendix A, Step 4(f): update z
	    # # # ###
	    
		# # # # p1 <- p*dmvt2(s,S.tilde[h,3:4],lc,Sigma,log=FALSE)
		# # # # p2 <- (1-p)*dmvt2(s,S.tilde[h,3:4],lc,Q,log=FALSE)
		# # # # p <- exp(log(p1)-log(p1+p2))	
		# # # # z <- rbinom(Ts,1,p)

		# # # ###
		# # # ### Appendix A, Step 4(g): update sigma2.alpha
		# # # ###
		
		# # # # r.tmp <- 1/(sum((alpha-mu.alpha)^2)/2+1/priors$r.sigma.alpha)
		# # # # q.tmp <- qW/2+priors$q.sigma.alpha
		# # # # sigma2.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		# # # # diag(Sigma.alpha) <- sigma2.alpha
		# # # # Sigma.alpha.inv <- solve(Sigma.alpha)
		# # # # A.inv.alpha <- solve(W.cross+Sigma.alpha.inv)
		

	    # # # #--------------------------------------------------------------------------
	  	# # # # Appendix A, Step 5: update sigma.mu
	    # # # #--------------------------------------------------------------------------
# # # # browser()

# # # priors$mu.sigma <- 6500/s.sd  # variance on lognormal prior for sigma.mu
# # # priors$tau <- 0.1

# # # sigma.mu.star <-  exp(rnorm(1,log(sigma.mu),tune$sigma.mu))
		# # # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,a,rho)
		# # # idx <- which(z==0)
	    # # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))					+dnorm(log(sigma.mu.star),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))						    +dnorm(log(sigma.mu),log(priors$mu.sigma),priors$tau,log=TRUE)		    
	    # # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
        	# # # sigma.mu <- sigma.mu.star
			# # # Q <- Q.star
			# # # keep$sigma.mu <- keep$sigma.mu+1
			# # # # keep.tmp$sigma.mu <- keep.tmp$sigma.mu+1

  	# # # } 



	    # # # # sigma.mu.star <- rnorm(1,sigma.mu,tune$sigma.mu)
	    # # # # if(sigma.mu.star>priors$sigma.mu.l & sigma.mu.star<priors$sigma.mu.u){
			# # # # # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,n.lc,Mix)
			# # # # Q.star <- get.Sigma(sigma2+sigma.mu.star^2,a,rho)
			# # # # idx <- which(z==0)
		    # # # # mh.star.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q.star,log=TRUE))
		    # # # # mh.0.sigma.mu <- sum(dmvt2(s[idx,],S.tilde[h[idx],3:4],lc[idx],Q,log=TRUE))
		    # # # # if(exp(mh.star.sigma.mu-mh.0.sigma.mu)>runif(1)){
	        	# # # # sigma.mu <- sigma.mu.star
				# # # # Q <- Q.star
				# # # # keep$sigma.mu <- keep$sigma.mu+1
	    	# # # # } 
	    # # # # }
		
				
		# # # #--------------------------------------------------------------------------
		# # # # Update spatial haul-out process model parameters: Appendix A, Step 6
		# # # #--------------------------------------------------------------------------
# # # # browser()
		# # # ###
		# # # ### Appendix A, Step 6(a): tabulate cluster membership 
		# # # ###

		# # # n <- table(h)  
		# # # m <- length(n)  # number of clusters
		# # # mu <- as.numeric(names(n))  # idx of occupied clusters

		# # # ###
		# # # ### Appendix A, Step 6(d): update 'unoccupied' mu
		# # # ###

		# # # p <- exp(U[-mu,]%*%gamma)
		# # # mu.tmp <- sample(S.tilde[-mu,1],J-m,replace=FALSE,prob=p)  # idx of unoccupied mu

		# # # ###
		# # # ### Appendix A, Step 6(c): update 'occupied' mu		   
		# # # ###
# # # # browser()			
		# # # mu.star <- sapply(mu,function(x) sample(S.tilde[,1],1,prob=  # proposals for mu	
			# # # dnorm(S.tilde[x,3],S.tilde[,3],tune$mu)*	dnorm(S.tilde[x,4],S.tilde[,4],tune$mu)))  
		# # # dup.idx <- which(!duplicated(mu.star))  # exclude duplicate proposals
		# # # mh <- sapply(dup.idx,function(x)  # accepted proposals
			# # # get.mh.mu(x,mu,mu.star,mu.tmp,S.tilde,h,s,lc,z,Sigma,Q,U,gamma)) 
		# # # keep$mu <- keep$mu+sum(mh)
		# # # keep.idx <- dup.idx[mh]
		# # # mu[keep.idx] <- mu.star[keep.idx]  # update mu

		# # # ###
	    # # # ### Appendix A, Step 6(f): update gamma
	    # # # ###

		# # # # Integral over S.tilde, occupied mu only
		# # # gamma.star <- matrix(rnorm(qU,gamma,tune$gamma),qU)
		# # # mh.star.gamma <- sum(dnorm(gamma.star,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # # sum(U[mu,]%*%gamma.star-log(sum(exp(U%*%gamma.star))))
 		 	# # # # sum(U[mu,]%*%gamma.star-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma.star))))
			# # # mh.0.gamma <- sum(dnorm(gamma,mu.gamma,priors$sigma.gamma,log=TRUE))+
			# # # sum(U[mu,]%*%gamma-log(sum(exp(U%*%gamma))))
			# # # # sum(U[mu,]%*%gamma-log(sum(exp(U[c(mu,mu.tmp),]%*%gamma))))
		# # # if(exp(mh.star.gamma-mh.0.gamma)>runif(1)){
    	    # # # gamma <- gamma.star
	        # # # # gamma.int <- gamma.int.star
	        # # # keep$gamma <- keep$gamma+1
    	# # # } 

		# # # ###
		# # # ### Appendix A, Step 6(b): update the stick-breaking process 
		# # # ###
		
			# # # # Appendix A, Step 6(b(i)): create index set I
			# # # I <- order(n,decreasing=TRUE)  # clusters ordered by membership
		
			# # # # Appendix A, Step 6(b(ii)): update eta
			# # # n.tmp <- c(n[I],rep(0,J-m-1))  # membership in decreasing order
			# # # eta <- c(rbeta(J-1,1+n.tmp,theta+Ts-cumsum(n.tmp)),1)  # stick-breaking weights

			# # # # Appendix A, Step 6(b(ii)): update pi
		    # # # pie <- eta*c(1,cumprod((1-eta[-J])))  # mixture component probabilities

		    # # # # Appendix A, Step 6(b(iii)): update theta
			# # # theta <- rgamma(1,priors$r.theta+J-1,priors$q.theta-sum(log(1-eta[-J])))  

			
		# # # ###
	    # # # ### Appendix A, Step 6(e): update h(t_s)
	    # # # ###

		# # # mu <- c(mu[I],mu.tmp)
		# # # h <- sapply(1:Ts,function(x) sample(mu,1,prob= 
			# # # exp(log(pie)+get.h(s[x,],S.tilde[mu,3:4],z[x],lc[x],Sigma,Q))))

		# # # ###
		# # # ###  Appendix A, Step 7: save samples 		   
		# # # ###
		
		# # # mu.save[,k] <- S.tilde[h,2]
		# # # theta.save[k] <- theta    
		# # # sigma.mu.save[k] <- sigma.mu
		# # # # sigma.alpha.save[k] <- sqrt(sigma2.alpha)
		# # # # alpha.save[k,] <- alpha
		# # # beta.save[k,] <- beta
		# # # gamma.save[k,] <- gamma
		# # # v.save[,k] <- v
		# # # z.save[,k] <- z
		# # # m.save[k] <- m
	# # # }
  	
  	# # # sigma.mu.save <- sigma.mu.save*s.sd
  	# # # t.mcmc.end <- Sys.time()

	# # # #####################################################################
	# # # ### Write output
	# # # #####################################################################
	  
	# # # keep$sigma.mu <- keep$sigma.mu/n.mcmc
	# # # keep$mu <- keep$mu/sum(m.save)
	# # # keep$gamma <- keep$gamma/n.mcmc
	# # # cat(paste("\n\ngamma acceptance rate:",round(keep$gamma,2))) 
	# # # cat(paste("\nmu acceptance rate:",round(keep$mu,2))) 
	# # # cat(paste("\nsigma.mu acceptance rate:",round(keep$sigma.mu,2))) 
	# # # cat(paste("\n\nEnd time:",Sys.time()))
	# # # cat(paste("\nTotal time elapsed:",round(difftime(Sys.time(),t.start,units="mins"),2)))
	# # # cat(paste("\nTime per MCMC iteration:",
		# # # round(difftime(t.mcmc.end,t.mcmc.start,units="secs")/n.mcmc,2)),"seconds")
	# # # cat(paste("\nTime per v update:",round(t.v.update/n.mcmc,2),"seconds"))
	# # # list(beta=beta.save,gamma=gamma.save,  # alpha=alpha.save,
		# # # mu=mu.save,theta=theta.save,m=m.save,z=z.save,v=v.save,
		# # # sigma.mu=sigma.mu.save,  # sigma.alpha=sigma.alpha.save,
		# # # keep=keep,n.mcmc=n.mcmc)
# # # }
