rm(list=ls())

# library(cluster)
library(msm)  # for simulating mu
# library(pscl)
library(MCMCpack)  # for rdirichlet(...)


###
### Simulate 2-dimensional haul out data using a Dirichlet process mixture 
###

# Define 2-dimensional rectangular support for mu, the true harbor seal locations
S.bar <- cbind(c(-10,-0.5,-0.5,-10,-10),c(0,0,20,20,0))  # complement of S, i.e., the land domain
S.tilde <- cbind(c(max(S.bar[,1]),0,0,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support 
	# of haul out process
S <- cbind(c(max(S.bar[,1]),10,10,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support of
	# movement process (marine and haul-out environments)

# Simulate cluster locations and assignments using stick-breaking process 
# See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
T <- 500  # number of locations to simulate
a0 <- 1.0  # concentration parameter
H <- 25  # maximum number of clusters for truncation approximation

mu.0 <- cbind(runif(H,min(S.tilde[,1]),max(S.tilde[,1])),
	runif(H,min(S.tilde[,2]),max(S.tilde[,2])))  # clusters randomly drawn from S.tilde
v <- c(rbeta(H-1,1,a0),1)  # stick-breaking weights
pie <- v*c(1,cumprod((1-v[-H])))  # probability mass
h.idx <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
h <- mu.0[h.idx,]  # latent clusters

# Simulate true locations
p <- 0.5  # probability of being hauled-out
z <- rbinom(T,1,p)  # haulout indicator variable: 1=hauled-out, 0=at-sea
table(z)
sigma.mu <- 2  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)
mu[z==1,] <- h[z==1,]
mu[z==0,] <- cbind(rtnorm(T-sum(z),h[z==0,1],sigma.mu,lower=S[1,1],upper=S[2,1]),
	rtnorm(T-sum(z),h[z==0,2],sigma.mu,lower=S[1,2],upper=S[3,2]))

# S.tilde <- S
# p <- 1
# z <- rep(1,T)
# mu <- h


# Observation process
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
points(mu,pch=19,cex=0.5) # All true locations
points(s,col=2,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6)) # Haul out locations


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulout.dp.mixture.mcmc.R")
start <- list(a0=a0,h=h,mu=mu,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r=1,q=0.25,sigma.l=0,sigma.u=5,sigma.mu.l=0,sigma.mu.u=5,
	alpha=1,beta=1)
# hist(rgamma(1000,2,0.5))
out1 <- haulout.dpmixture.mcmc(s,S.tilde,S,priors=priors,
  tune=list(z=0.5,mu.0=0.15,sigma=0.05,sigma.mu=0.25),start=start,n.mcmc=2000)

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
points(s,pch=19,cex=0.2,col=3)
points(h,pch=1,cex=1,col=rgb(1,0,0,1))

pt.idx <- 4
points(mod$mu[pt.idx,1,idx],mod$mu[pt.idx,2,idx],pch=19,cex=0.2,col=rgb(0,1,0,0.25))
points(mu[pt.idx,1],mu[pt.idx,2],pch=19,cex=0.57)
points(s[pt.idx,1],s[pt.idx,2],pch=19,col=2)
table(mod$z[pt.idx,idx])
z[pt.idx]

# Concentration parameter
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2) 
mean(mod$a0[idx])*log(T)

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)
mean(mod$sigma[idx])

# Dispersion about haul-out for at-sea locations
hist(mod$sigma.mu[idx],breaks=100);abline(v=sigma.mu,col=2,lty=2)
mean(mod$sigma.mu[idx])

# Haul-out probability
hist(mod$p[idx],breaks=100);abline(v=p,col=2,lty=2)
abline(v=sum(z)/T,col=3,lty=2)

# Haul-out indicator variable
z.hat <- apply(mod$z[,idx],1,sum)/(length(idx))
boxplot(z.hat~z,ylim=c(0,1))
plot(s[,1],z.hat,ylim=c(0,1))

# Modeled number of clusters
# n.cls <- apply(mod$h[,,idx],3,function(x) nrow(unique(x)))
plot(mod$n.cls,type="l")
abline(h=nrow(unique(h)),col=2,lty=2)  # true number of clusters  
barplot(table(mod$n.cls))





plot(apply(mod$z[,idx],2,max),type="l")
cl.ranks <- apply(mod$z[,idx],2,dense_rank)
plot(c(mod$z[,idx])[c(cl.ranks)==8],type="l")

hist(c(mod$z[,idx])[c(test)==1],breaks=500,xlim=range(mod$z[idx]),ylim=c(0,5),prob=TRUE)  
hist(c(mod$z[,idx])[c(test)==5],col=5,breaks=500,add=TRUE,prob=TRUE)  

pt.idx <- 94
plot(mod$z[pt.idx,idx],type="l");abline(h=z[pt.idx],col="red",lty=2)
hist(mod$z[,idx],breaks=5000,xlim=c(range(mod$z[pt.idx,])+c(-10,10)),prob=TRUE)
hist(mod$z[pt.idx,idx],breaks=50,col="red",add=TRUE,border="red",prob=TRUE);abline(v=z[pt.idx],col="red",lty=2)
points(y[pt.idx],-0.010,pch=19)

which.max(apply(mod$z,1,var))
abline(v=29.75)
which.min(sapply(y,function(x) dist(c(x,29.75))))
