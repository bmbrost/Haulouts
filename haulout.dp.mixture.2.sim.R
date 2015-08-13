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

S.bar <- cbind(c(-10,-0.5,-0.5,-10,-10),c(0,0,20,20,0))  # complement of S, i.e., the land domain
S.tilde <- cbind(c(max(S.bar[,1]),0,0,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support 
	# of haul out process
S <- cbind(c(max(S.bar[,1]),10,10,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support of
	# movement process (marine and haul-out environments)

###
### Simulate cluster locations and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
###

T <- 500  # number of locations to simulate
a0 <- 1.0  # concentration parameter
H <- 25  # maximum number of clusters for truncation approximation
mu.0 <- cbind(runif(H,min(S.tilde[,1]),max(S.tilde[,1])),
	runif(H,min(S.tilde[,2]),max(S.tilde[,2])))  # clusters randomly drawn from S.tilde
v <- c(rbeta(H-1,1,a0),1)  # stick-breaking weights
pie <- v*c(1,cumprod((1-v[-H])))  # probability mass
h.idx <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
h <- mu.0[h.idx,]  # latent clusters

###
### Simulate wet/dry status for T telemetry locations (s) and T SEA records (y)
###

time <- c(0,cumsum(rgamma(T*2-1,shape=1.1,scale=2)))  # time covariate
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
s.idx <- sort(sample(1:(T*2),T))  # subset times between telemetry locs (s) and SEA data (y)
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
sigma.mu <- 2  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)
mu[z==1,] <- h[z==1,]
mu[z==0,] <- cbind(rtnorm(T-sum(z),h[z==0,1],sigma.mu,lower=S[1,1],upper=S[2,1]),
	rtnorm(T-sum(z),h[z==0,2],sigma.mu,lower=S[1,2],upper=S[3,2]))

# Simulate observed locations
sigma <- 0.5 # Observation error
s <- mu
s <- s+rnorm(T*2,0,sigma) # Add error to true locations

# Plot support
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)

# Plot true and observed locations
points(mu,pch=19,cex=0.25,col=z+1) # All true locations
points(s,col=z+1,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6),cex=0.5) # Haul out locations


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulout.dp.mixture.2.mcmc.R")
start <- list(a0=a0,h=h,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie,beta=beta)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r=4,q=2,sigma.l=0,sigma.u=5,sigma.mu.l=0,sigma.mu.u=5,
	alpha=1,beta=1,sigma.beta=1)
tune <- list(sigma=0.05,sigma.mu=0.25,a0=0.25)
# hist(rgamma(1000,4,2))
# hist(rgamma(1000,5,2.5))
out1 <- haulout.dpmixture.2.mcmc(s,y,X[s.idx,],X[-s.idx,],W[s.idx,],W[-s.idx,],
	S,S.tilde,sigma.alpha=2,priors=priors,tune=tune,start=start,n.mcmc=2000)

mod <- out1 
# idx <- 1:100
idx <- 1:1000
idx <- 1:2000
idx <- 1:10000

# True clusters
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)
points(mod$h[,1,idx],mod$h[,2,idx],pch=19,cex=0.5,col=rgb(0,0,0,0.0025))
points(h,pch=1,cex=1,col=rgb(1,0,0,1))
points(s,pch=19,cex=0.2,col=3)

pt.idx <- 101
points(mod$h[pt.idx,1,idx],mod$h[pt.idx,2,idx],pch=19,cex=0.2,col=rgb(1,1,0,0.25))
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

# Concentration parameter
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2) 
mean(mod$a0[idx])*log(T)

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)
mean(mod$sigma[idx])

# Dispersion about haul-out for at-sea locations
hist(mod$sigma.mu[idx],breaks=100);abline(v=sigma.mu,col=2,lty=2)
mean(mod$sigma.mu[idx])

# Inference on z: latent haul-out indicator variable for telemetry locations
z.hat <- apply(mod$z[,idx],1,sum)/(length(idx))
boxplot(z.hat~z,ylim=c(0,1))
plot(s[,1],z.hat,ylim=c(0,1),col=z+1)

u <- apply(out1$beta[idx,],1,function(x) X[s.idx,]%*%x)+
	apply(out1$alpha[idx,],1,function(x) W[s.idx,]%*%x)
u.inv <- matrix(pnorm(u),,T,byrow=TRUE)
u.inv.mean <- apply(u.inv,2,mean)
plot(u.inv.mean,z.hat)
u.inv.quant <- t(apply(u.inv,2,quantile,c(0.025,0.975)))
plot(u.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,u.inv.quant[,1],1:T,u.inv.quant[,2],col=rgb(0,0,0,0.15))
abline(h=0.5,col=2,lty=2)
points(z,col=3,pch=19,cex=0.5)
points(z.hat,pch=19,cex=0.25)

# Inference for y
u.tilde.inv <- matrix(pnorm(out1$u[,idx]),,T,byrow=TRUE)
u.tilde.inv.mean <- apply(u.tilde.inv,2,mean)
u.tilde.inv.quant <- t(apply(u.tilde.inv,2,quantile,c(0.025,0.975)))
plot(u.tilde.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,u.tilde.inv.quant[,1],1:T,u.tilde.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y,col=4,pch=19,cex=0.5)

# Modeled number of clusters
# n.cls <- apply(mod$h[,,idx],3,function(x) nrow(unique(x)))
plot(mod$n.cls,type="l")
abline(h=nrow(unique(h)),col=2,lty=2)  # true number of clusters  
barplot(table(mod$n.cls)) 


