##################################################################################
### Simulate 2-dimensional haul out data using a Dirichlet process mixture 
##################################################################################

rm(list=ls())

# library(cluster)
library(msm)  # for simulating mu
# library(pscl)
library(MCMCpack)  # for rdirichlet(...)
library(splines)

###
### Define 2-dimensional rectangular support for mu, the true harbor seal locations
###

S.bar <- cbind(c(-10,-0.5,-0.5,-10,-10),c(0,0,10,10,0))  # complement of S, i.e., the land domain
n.haulout <- 1000
S.tilde <- cbind(max(S.bar[,1]),  # support of haul out process
	seq(min(S.bar[,2]),max(S.bar[,2]),length.out=n.haulout)) 
S <- cbind(c(max(S.bar[,1]),10,10,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support of
	# movement process (marine and haul-out environments)

###
### Simulate cluster locations and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
###

T <- 500  # number of locations to simulate
n <- 500  # number of wet/dry observations to simulate 
theta <- 2.0  # concentration parameter
H <- 25  # maximum number of clusters for truncation approximation

# Resource selection strength for haul-out sites
U <- cbind(1,S.tilde[,2],sin(S.tilde[,2]))  # design matrix for haul-out site RSF
U[,-1] <- scale(U[,-1])
gamma <- c(-1,0.025,0.5)  # Haul-out site selection coefficients

U <- cbind(1,sin(5*S.tilde[,2]))  # design matrix for haul-out site RSF
U[,-1] <- scale(U[,-1])
gamma <- c(0,3)  # Haul-out site selection coefficients

# U <- cbind(1,rep(c(0,1),50))
# gamma <- c(0,2)  # Haul-out site selection coefficients

qU <- ncol(U)
selection <- exp(U%*%gamma)
plot(S.tilde[,2],selection,type="l")

# idx <- order(selection,decreasing=TRUE)
# plot(selection[idx],type="l")
# selection[-c(1:4)] <- 0

mu.0 <- sample(1:n.haulout,H,replace=TRUE,prob=selection)
mu.0 <- S.tilde[mu.0,]
length(unique(mu.0[,2]))==H
eta <- c(rbeta(H-1,1,theta),1)  # stick-breaking weights
pie <- eta*c(1,cumprod((1-eta[-H])))  # probability mass
h.match <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
h <- mu.0[h.match,]  # latent clusters

plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1])),ylim=range(S.bar[,2]),
	pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
points(S.tilde,pch=19,cex=0.25)
points(mu.0)
points(h,col=2,pch=19,cex=0.5)


###
### Simulate wet/dry status for T telemetry locations (s) and T SEA records (y)
###

time <- c(0,cumsum(rgamma(T+n-1,shape=1.1,scale=2)))  # time covariate
hr <- ((time)-24*floor(time/24))
day <- ceiling(time/24)
X <- cbind(1,day,hr)  # design matrix for telemetry locations and SEA data
qX <- ncol(X)

# Center and scale design matrix
X[,-1] <- scale(X[,-1])

beta <- c(0.75,1.75,1.0)  # Coefficients on X

trend <- 2*sin(0.01*time)  # non-linear pattern
plot(time,trend,type="l")

# Calculate probability of record (telemetry or SEA) being hauled-out
p <- pnorm(X%*%beta+trend)  # probability of being hauled-out
hist(p);summary(p)

# Simulate wet/dry status for telemetry locations (s) and SEA data (y)
s.idx <- sort(sample(1:(T+n),T))  # subset times between telemetry locs (s) and SEA data (y)
z <- rbinom(T,1,p[s.idx])  # Haul-out indicator variable for s: 1=hauled-out, 0=at-sea
table(z)
y <- rbinom(T,1,p[-s.idx])  # Haul-out indicator variable for y: 1=hauled-out, 0=at-sea
table(y)

# Inspect wet/dry status, along with covariates and trends
plot(c(time[s.idx],time[-s.idx]),c(z,y),pch="|",cex=0.5,col=c(rep(1,T),rep(2,T)),
	ylim=range(c(trend,y,X%*%beta)))
lines(time,X%*%beta,col=3)
lines(time,trend,col=4)
lines(time,X%*%beta+trend,col=5)

###
### B-splines basis expansion
###

int <- 50  # interval between knots
knots <- seq(0,max(time),by=int)
W <- bs(time,knots=knots,degree=3,intercept=FALSE)  # cubic spline
matplot(W,type="l")
qW <- ncol(W)

###
### Simulate true and observed locations
###

# Simulate true locations
sigma.mu <- 1  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)
mu[z==1,] <- h[z==1,]
mu[z==0,] <- cbind(rtnorm(T-sum(z),h[z==0,1],sigma.mu,lower=S[1,1],upper=S[2,1]),
	rtnorm(T-sum(z),h[z==0,2],sigma.mu,lower=S[1,2],upper=S[3,2]))

# Simulate observed locations
sigma <- 0.5 # Observation error
s <- mu
s <- s+rnorm(T*2,0,sigma) # Athetadd error to true locations

# Plot support
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
points(S.tilde,pch=19,cex=0.25)

# Plot true and observed locations
points(mu,pch=19,cex=0.25,col=z+1) # All true locations
points(s,col=z+1,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6),cex=0.5) # Haul out locations


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulout.dp.mixture.3.mcmc.R")
start <- list(theta=theta,h=h,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie,beta=beta,gamma=gamma)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r=4,q=2,sigma.l=0,sigma.u=5,sigma.mu.l=0,sigma.mu.u=5,
	sigma.beta=5,mu.gamma=rep(0,qU),sigma.gamma=5)
tune <- list(sigma=0.05,sigma.mu=0.2,mu.0=20,gamma=3.0)
# hist(rgamma(1000,4,2))
# hist(rgamma(1000,5,2.5))
out1 <- haulout.dpmixture.3.mcmc(s,y,X[s.idx,],X[-s.idx,],W[s.idx,],W[-s.idx,],
	S,S.tilde,U,sigma.alpha=2,priors=priors,tune=tune,start=start,n.mcmc=20000)

mod <- out1
idx <- 1:2000
idx <- 1:10000
idx <- 1:20000

# True clusters
h.tab <- table(mod$h[,2,idx])
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
# polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)
points(rep(-0.5,length(h.tab)),as.numeric(names(h.tab)),pch=19,cex=h.tab/max(h.tab),col=5)
points(S.tilde+cbind(rep(-1.5,n.haulout),0),cex=selection/max(selection),pch=19)
# points(mod$h[,1,idx],mod$h[,2,idx],pch=19,cex=0.5,col=rgb(1,0,0,0.001))
points(h,pch=1,cex=1,col=rgb(1,0,0,1))
points(s,pch=19,cex=0.2,col=3)

plot(rep(-0.5,length(h.tab)),as.numeric(names(h.tab)),pch=19,cex=h.tab/max(h.tab),
	ylim=range(S.tilde),asp=TRUE)
points(S.tilde+cbind(rep(-1.5,n.haulout),0),cex=selection/max(selection),pch=19)
# plot(jitter(mod$h[,1,idx]),jitter(mod$h[,2,idx]),
	# pch=19,cex=0.5,col=rgb(0,0,0,0.0025),asp=TRUE,ylim=range(S.tilde))
points(h,pch=1,cex=1,col=rgb(1,0,0,1))
points(s,pch=19,cex=0.2,col=3)

idx <- 850
h.tab[which(names(h.tab)==S.tilde[idx,2])]
points(S.tilde[idx,1],S.tilde[idx,2])
junk <- mod$gamma[which(apply(mod$h,3,function(x) sum(x[,2]==S.tilde[idx,2]))>0),]
plot(S.tilde[,2],exp(U%*%gamma),type="l")

max(exp(U%*%gamma))/min(exp(U%*%gamma))
apply(junk,1,function(x) max(exp(U%*%x))/min(exp(U%*%x)))
plot(S.tilde[,2],exp(U%*%junk[12,]),type="l")

junk <- apply(mod$gamma,1,function(x) max(exp(U%*%x))/min(exp(U%*%x)))
summary(junk)
which(junk<1.1)

hist(junk[junk<50],breaks=100)

which(s[,1] < (-1.5))
pt.idx <- 37
points(mod$h[pt.idx,1,idx],mod$h[pt.idx,2,idx],pch=19,cex=0.2,col=rgb(1,1,0,0.5))
points(mu[pt.idx,1],mu[pt.idx,2],pch=19,cex=0.75)
points(s[pt.idx,1],s[pt.idx,2],pch=19,col=2)
table(mod$z[pt.idx,idx])
z[pt.idx]

# Haul-out probability covariates
matplot(mod$beta[idx,],type="l"); abline(h=beta,col=1:3,lty=2) 
beta.hat <- apply(out1$beta[idx,],2,mean)
beta.quant <- t(apply(out1$beta[idx,],2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)

# Haul-out selection covariates
matplot(mod$gamma[idx,-1],type="l"); abline(h=gamma[-1],col=1:3,lty=2) 
gamma.hat <- apply(out1$gamma[idx,],2,mean)
gamma.quant <- t(apply(out1$gamma[idx,],2,quantile,c(0.025,0.975)))
plot(gamma.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(gamma.quant)))
abline(h=0,col=2,lty=2)
segments(1:qU,gamma.quant[,1],1:qU,gamma.quant[,2],col="lightgrey")
points(gamma.hat,pch=19,col=rgb(0,0,0,0.25))
points(gamma,pch=19)

# Inference on alpha
alpha.hat <- apply(out1$alpha[idx,],2,mean)
alpha.quant <- t(apply(out1$alpha[idx,],2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:qW,alpha.quant[,1],1:qW,alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
# points(alpha,pch=19,col=3)

par(mfrow=c(2,1))
plot(time,trend,type="l")
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
lines(alpha.hat,col=rgb(0,0,0,0.25))
abline(h=0,col=2,lty=2)

# Concentration parameter (theta) and number of clusters
hist(mod$theta[idx],breaks=100);abline(v=theta,col=2,lty=2) 
mean(mod$theta[idx])*log(T)

# n.cls <- apply(mod$h[,,idx],3,function(x) nrow(unique(x)))
plot(mod$n.cls,type="l");abline(h=nrow(unique(h)),col=2,lty=2)  # true number of clusters  
barplot(table(mod$n.cls)) 

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)
mean(mod$sigma[idx])

# Dispersion about haul-out for at-sea locations
hist(mod$sigma.mu[idx],breaks=100);abline(v=sigma.mu,col=2,lty=2)
mean(mod$sigma.mu[idx])

# Inference on z: latent haul-out indicator variable for telemetry locations
# Note: estimation of z (z.hat) is based on covariates and location of telemetry 
# observations, whereas v is based on covariates alone.
z.hat <- apply(mod$z[,idx],1,sum)/(length(idx))
boxplot(z.hat~z,ylim=c(0,1))
plot(s[,1],z.hat,ylim=c(0,1),col=z+1);abline(v=S.tilde[1,1],col=3,lty=2)

v <- apply(mod$beta[idx,],1,function(x) X[s.idx,]%*%x)+
	apply(mod$alpha[idx,],1,function(x) W[s.idx,]%*%x)
v.inv <- matrix(pnorm(v),,T,byrow=TRUE)
v.inv.mean <- apply(v.inv,2,mean)
plot(v.inv.mean,z.hat)
v.inv.quant <- t(apply(v.inv,2,quantile,c(0.025,0.975)))
plot(v.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,v.inv.quant[,1],1:T,v.inv.quant[,2],col=rgb(0,0,0,0.15))
abline(h=0.5,col=2,lty=2)
points(z,col=3,pch=19,cex=0.5)
points(z.hat,pch=19,cex=0.25,col=1)

# Inference for y
v.tilde.inv <- matrix(pnorm(mod$v[-(1:T),idx]),,n,byrow=TRUE)
v.tilde.inv.mean <- apply(v.tilde.inv,2,mean)
v.tilde.inv.quant <- t(apply(v.tilde.inv,2,quantile,c(0.025,0.975)))
plot(v.tilde.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,v.tilde.inv.quant[,1],1:T,v.tilde.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y,col=4,pch=19,cex=0.5)
