rm(list=ls())

# library(cluster)
library(msm)  # for simulating mu
# library(pscl)
library(MCMCpack)  # for rdirichlet(...)


###
### Simulate 2-dimensional haul out data using a Dirichlet process mixture 
###

# Define 2-dimensional rectangular support for mu, the true harbor seal locations
S.bar <- cbind(c(-10,-2,-2,-10,-10),c(0,0,20,20,0))  # complement of S, i.e., the land domain
S.tilde <- cbind(c(max(S.bar[,1]),1,1,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support 
	# of haul out process
S <- cbind(c(max(S.bar[,1]),10,10,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2])  # support of
	# movement process (marine and haul-out environments)

# Simulate cluster locations and assignments using stick-breaking process 
# See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
T <- 500  # number of locations to simulate
a0 <- 1.5  # concentration parameter
H <- 50  # maximum number of clusters for truncation approximation

theta <- cbind(runif(H,min(S.tilde[,1]),max(S.tilde[,1])),
	runif(H,min(S.tilde[,2]),max(S.tilde[,2])))  # clusters randomly drawn from S.tilde
v <- c(rbeta(H-1,1,a0),1)  # stick-breaking weights
pie <- v*c(1,cumprod((1-v[-H])))  # probability mass
h <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
h <- theta[h,]  # latent clusters

# Simulate true locations
p <- 0.5  # probability of being hauled-out
z <- rbinom(T,1,p)  # haulout indicator variable
sigma.mu <- 2  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)
mu[z==1,] <- h[z==1,]
mu[z==0,] <- cbind(rtnorm(T-sum(z),h[z==0,1],sigma.mu,lower=S[1,1],upper=S[2,1]),
	rtnorm(T-sum(z),h[z==0,2],sigma.mu,lower=S[1,2],upper=S[3,2]))

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
points(mu) # All true locations
points(s,col=2,pch=3) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6)) # Haul out locations


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulout.dp.mixture.mcmc.R")
start <- list(a0=a0,h=h,mu=mu,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie)  # rdirichlet(1,rep(1/H,H))) 
out1 <- haulout.dpmixture.mcmc(s,S.tilde,S,
  priors=list(H=H,r=20,q=10,sigma.l=0,sigma.u=5,sigma.mu.l=0,sigma.mu.u=5),
  tune=list(z=0.5,sigma=0.01,sigma.mu=0.01),start=start,n.mcmc=1000)

mod <- out1
idx <- 1:100
idx <- 1:1000
idx <- 1:2500
idx <- 1:10000

# True clusters
b <- 3*c(-sigma,sigma) # Plot buffer for errors
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)
points(mod$h[,1,idx],mod$h[,2,idx],pch=19,cex=0.5,col=rgb(0,0,0,0.025))
points(s,pch=19,cex=0.2,col=3)
points(h,pch=19,cex=0.5,col=rgb(1,0,0,1))

points(mod$mu[,1,idx],mod$mu[,2,idx],pch=19,cex=0.2,col=rgb(0,0,1,0.025))


cl <- kmeans(apply(mod$z,2,I),11)
points(cl$centers,col=4,pch=19)


# Concentration parameter
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2) 
mean(mod$a0[idx])*log(n)

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)

# Modeled number of clusters
nclust <- apply(mod$z[,,idx],c(3),function(x) nrow(unique(x)))
head(nclust)
plot(nclust,type="l")
abline(h=nrow(unique(z)),col=2,lty=2)  # true number of clusters  

barplot(table(nclust))




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
