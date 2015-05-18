rm(list=ls())

library(msm)
library(pscl)

###
### Simulate 2-dimensional haul out data 
###

# Define 2-dimensional rectangular support for mu
S.bar <- cbind(c(-10,-2,-2,-10,-10),c(0,0,20,20,0)) # Complement of S, i.e., the land domain
S.tilde <- cbind(c(max(S.bar[,1]),2,2,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2]) # Support of haul out process
  c(max(S.bar),2) 
S <- cbind(c(max(S.bar[,1]),10,10,max(S.bar[,1]),max(S.bar[,1])),S.bar[,2]) # Support of movement process (marine and haul-out environments)

# Biological process
T <- 100 # Number of locations
p <- 0.5 # Probability of being on the haul-out
z <- rbinom(T,1,p)
mu.0 <- c(runif(1,S.tilde[1,1],S.tilde[2,1]),runif(1,S.tilde[1,2],S.tilde[3,2])) # Homerange center
sigma.mu <- 5

mu <- matrix(0,T,2)
mu[z==1,] <- cbind(runif(sum(z),S.tilde[1,1],S.tilde[2,1]),runif(sum(z),0,20))
mu[z==0,] <- cbind(rtnorm(T-sum(z),mu.0[1],sigma.mu,lower=S[1,1],upper=S[2,1]),
                   rtnorm(T-sum(z),mu.0[2],sigma.mu,lower=S[1,2],upper=S[3,2]))

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
points(mu.0[1],mu.0[2],pch=17)
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6)) # Haul out locations
points(s,col=2,pch=3) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu


###
### Fit model
###

source("/Users/brost/Documents/git/Haulouts/haulout.homerange.mcmc.R")
priors <- list(alpha=1,beta=1,q=3,r=2,a=0,b=2)
tune <- list(mu=0,sigma=0.1,sigma.mu=0.5)
start <- list(mu=mu,sigma=sigma,mu.0=mu.0,sigma.mu=sigma.mu)
out1 <- haulout.homerange.mcmc(s,S.tilde,S,priors,tune,start,n.mcmc=10000)

mod <- out1
plot(mod$p,type="l");abline(h=p,col=2)
abline(h=sum(z)/T,col=2)
abline(h=mean(mod$p),lty=2,col=3)
plot(mod$sigma,type="l");abline(h=sigma,col=2)
table(apply(mod$z,1,sum)>out1$n.mcmc/2)

sum(z)

idx <- mu[,1]>S.tilde[1,1]&mu[,1]<S.tilde[2,1]&
  mu[,2]>S.tilde[1,2]&mu[,2]<S.tilde[3,2] #mu located within intersection(S,S.tilde)
sum(idx)

# Posterior of mu[t] and z[t]
mu.hat <- apply(mod$mu,c(1,2),mean)
idx <- 39
z[idx]
table(mod$z[idx,])
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)
points(t(mod$mu[idx,,]),col=rgb(0,0,0,0.1),pch=19)
points(mu.hat[idx,1],mu.hat[idx,2],col=3,pch=19)
points(mu[idx,1],mu[idx,2],col=2,pch=19)
points(s[idx,1],s[idx,2],col=4,pch=19)


# Posterior of mu.0
plot(0,0,xlim=c(min(S.bar[,1]),max(S[,1]))+b,ylim=range(S.bar[,2])+b,pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=S.bar[,1],y=S.bar[,2],col="gray45")
polygon(x=S[,1],y=S[,2],col="gray85")
polygon(x=S.tilde[,1],y=S.tilde[,2],angle=45,density=5)

points(mod$mu.0,col=rgb(0,0,0,0.1),pch=19)
mu.0.hat <- apply(mod$mu.0,2,mean)
points(mu.0.hat[1],mu.0.hat[2],col=3,pch=19)


str(mod)




###
### Explore IG parameterization
###

# Inverse gamma using rigamma {pscl}
library(pscl)
igammastrt <- function(m,v){
  alpha <-m^2/v+2
  beta <- m^3/v+m
  list(alpha=alpha,beta=beta)
}

igmn <- sigma^2
igvar <- 0.5
igammastrt(igmn,igvar)
hist(rigamma(1000,igammastrt(igmn,igvar)$alpha,igammastrt(igmn,igvar)$beta)) # arguments: n,alpha,beta or n,q,r
mean(rigamma(10000,igammastrt(igmn,igvar)$alpha,igammastrt(igmn,igvar)$beta))
var(rigamma(10000,igammastrt(igmn,igvar)$alpha,igammastrt(igmn,igvar)$beta))

hist(sqrt(rigamma(1000,40,10))) # arguments: n, q, r
abline(v=sigma,col=2)


# Inverse gamma using rgamma {base}
invgammastrt <- function(igmn,igvar){
  q <- 2+(igmn^2)/igvar
  r <- 1/(igmn*(q-1))
  list(r=r,q=q)
}

invgammastrt(0.1,sqrt(0.25))
hist(1/rgamma(1000,2.04,,9.62)) # mean=0.5, variance=1
var(1/rgamma(100000,2.02,,9.80))

hist(1/rgamma(1000,2.25,,1.6)) # mean=0.5, variance=1
hist(1/rgamma(1000,2.04,,9.62)) # mean=0.1, variance=0.25

