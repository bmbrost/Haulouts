library(msm)
library(pscl)

###
### Simulate 1-dimensional haul out data 
###

# Define 1-dimensional support for mu
S.bar <- c(-10,-2) # Complement of S, i.e., the land domain
S.tilde <- c(max(S.bar),2) # Support of haul out process
S <- c(max(S.bar),10) # Support of movement process (marine and haul-out environments)

# Biological process
n <- 100 # Number of locations
p <- 0.5 # Probability of being on the haul-out
z <- rbinom(n,1,p)
mu <- cbind(0,runif(n,0,1))
mu[z==1,1] <- runif(sum(z),S.tilde[1],S.tilde[2])
mu[z==0,1] <- runif(n-sum(z),S[1],S[2])

# Observation process
sigma <- 1 # Observation error
s <- mu
s[,1] <- s[,1]+rnorm(n,0,sigma) # Add error to true locations

# Plot support
plot(0,0,xlim=c(min(S.bar),max(S)),ylim=c(0,1),pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=c(S.bar,rev(S.bar),S.bar[1]),y=c(0,0,1,1,1),col="gray45")
polygon(x=c(S,rev(S),S[1]),y=c(0,0,1,1,1),col="gray85")
polygon(x=c(S.tilde,rev(S.tilde),S.tilde[1]),y=c(0,0,1,1,1),angle=45,density=5)

# Plot true and observed locations
points(mu) # All true locations
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6)) # Haul out locations
points(s,col=2,pch=3) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu


###
### Fit model
###

source("haulout.1d.mcmc.R")
out1 <- haulout.1d.mcmc(s[,1],mu[,1],S.tilde,S,priors=list(alpha=1,beta=1,q=3,r=2,a=0,b=2),
  tune=list(mu=0.25,sigma=0.25),n.mcmc=5000)

plot(out1$p,type="l");
abline(h=sum(z)/n,col=2)
abline(h=mean(out1$p),lty=2,col=3)
plot(out1$sigma,type="l");abline(h=sigma,col=2)
table(apply(out1$z,1,sum)>out1$n.mcmc/2)
sum(z)

# Posterior of mu[t] and z[t]
idx <- 90
z[idx]
table(out1$z[idx,])
hist(out1$mu[idx,],col="gray80",breaks=seq(S[1],S[2],0.5))
abline(v=mean(out1$mu[idx,]),col=3);abline(v=mu[idx,1],col=2);abline(v=s[idx,1],lty=2)
hist(rtnorm(out1$n.mcmc,mu[idx,1],sigma,lower=S[1],upper=S[2]),add=TRUE,density=5,angle=45,breaks=seq(S[1],S[2],0.5))


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

