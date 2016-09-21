rm(list=ls())

library(sp)
library(rgdal)
library(maptools)
library(raster)
library(mvtnorm)
library(rgeos)
library(gdistance)
library(DPpackage)
library(chron)

# See R/hualouts.script.R for data prep and creation of sim.data
setwd("~/Documents/research/harbor_seals/projects/central_places/")
load('~/Documents/research/harbor_seals/projects/central_places/R/sim/sim.Rdata')
ls()

get.s <- function(mu,lc,sigma,a,rho,nu){ #Simulate observed locations (s[t])
	library(mvtnorm)
	s <- mu
	T <- nrow(s)
	lc.tmp <- sort(unique(lc))
	Sigma <- lapply(lc.tmp, function(x) 
		sigma[x]^2*matrix(c(1,sqrt(a[x])*rho[x],sqrt(a[x])*rho[x],a[x]),2))
	tdist.idx <- sample(1:2,T,replace=TRUE) # Index for mixture t-distribution
	K <- matrix(c(-1,0,0,1),2)
	# browser()
	for(i in 1:T){
		lc.tmp <- lc[i]
		if(tdist.idx[i]==1){	
			# s[i,] <- mu[i,]+rmvt(1,sigma=Sigma[[lc.tmp]],df=nu[lc.tmp])
			s[i,] <- mu[i,]+rmvnorm(1,sigma=Sigma[[lc.tmp]])
		}
		if(tdist.idx[i]==2){
			# s[i,] <- mu[i,]+rmvt(1,sigma=K%*%Sigma[[lc.tmp]]%*%t(K),df=nu[lc.tmp])
			s[i,] <- mu[i,]+rmvnorm(1,sigma=K%*%Sigma[[lc.tmp]]%*%t(K))
		}
	}
	s
}


##################################################################################
### Estimated parameters from PV95KOD09
##################################################################################

theta <- as.numeric(read.csv("R/sim/params/theta.hat.csv"))
sigma.mu <- as.numeric(read.csv("R/sim/params/sigma.mu.hat.csv"))
beta <- as.numeric(unlist(read.csv("R/sim/params/beta.hat.csv")))
alpha <- as.numeric(unlist(read.csv("R/sim/params/alpha.hat.csv")))
sigma <- as.numeric(unlist(read.csv("R/sim/params/sigma.hat.csv")))[c(1,4,6)]
a <- as.numeric(unlist(read.csv("R/sim/params/a.hat.csv")))[c(1,4,6)]
rho <- as.numeric(unlist(read.csv("R/sim/params/rho.hat.csv")))[c(1,4,6)]
sigma.alpha <- as.numeric(read.csv("R/sim/params/sigma.alpha.hat.csv"))


##################################################################################
### Simulate haul-out sites
##################################################################################

J <- 50  # Upper threshold for truncation approximation for DP
J <- 3  

# Stick-breaking process
eta <- c(rbeta(J-1,1,theta),1)  # stick-breaking weights
pi <- eta*c(1,cumprod((1-eta[-J])))  # probability mass
plot(pi,type="b")

mu <- sample(S.tilde.idx,J,prob=f.hat)  # clusters drawn from [mu.0|gamma,S.tilde]
plot(S.tilde.xy,pch="")
plot(S.tilde.tmp,add=TRUE,maxpixels=100^10)
points(xyFromCell(S.tilde,mu),col=1,pch=19,cex=0.25)  # possible haul-out sites

# Simulate realized haul-out sites
h <- sample(mu,Ts,replace=TRUE,prob=pi)  # latent cluster assignments;
	# note this references cell in S.tilde, not record in mu as in Appendix A
n <- table(h)  # tabulate cluster membership
m <- length(n)  # number of clusters

# Plot realized haul-out locations
points(xyFromCell(S.tilde,as.numeric(names(n))),pch=19,cex=n/max(n)+0.25,col=2)


##################################################################################
### Simulate wet/dry status
##################################################################################

p <- pnorm(X.tmp%*%beta+W%*%alpha)  # probability of being hauled-out
hist(p)
y <- rbinom(Ts,1,p)  # haulout indicator variable: 1=hauled-out, 0=at-sea
table(y)


##################################################################################
### Simulate true animal locations
##################################################################################

mu.t <- matrix(0,Ts,2)  # x,y coordinates of true locations for all s

# Sample true locations for hauled-out observations ; i.e., true location=haul-out site
mu.t[y==1,] <- xyFromCell(S.tilde,h[y==1])

# Sample true locations for at-sea observations 
mu.buffer <- gBuffer(mu.start,width=30000)
S.mask <- mask(S,mu.buffer)
plot(S.mask)
idx <- which(values(S.mask)==1)
S.xy <- xyFromCell(S,idx)
mu.t[y==0,1] <- apply(xyFromCell(S.tilde,h[y==0]),1,function(x) sample(idx,1,prob=
	dnorm(x[1],S.xy[,1],sigma.mu)*dnorm(x[2],S.xy[,2],sigma.mu))) 
mu.t[y==0,] <- xyFromCell(S,mu.t[y==0,1])

# plot(xyFromCell(S.tilde,as.numeric(names(n))),pch=19,cex=n/max(n)+0.25,col=2)
plot(mu.t,pch="")
plot(S,add=TRUE)
points(mu.t)
points(xyFromCell(S.tilde,as.numeric(names(n))),pch=19,cex=n/max(n)+0.25,col=2)

# Check true location generation
# mu.tmp <- xyFromCell(S.tilde,h)
# S.test <- S
# S.test[idx] <- dnorm(mu.tmp[10,1],S.xy[,1],sigma.mu)*dnorm(mu.tmp[10,2],S.xy[,2],sigma.mu)
# plot(mu.t,pch="",cex=0.25)
# plot(S.tilde,add=TRUE)
# plot(S.test,add=TRUE)
# points(mu.tmp[10,1],mu.tmp[10,2])
# points(mu.t[10,1],mu.t[10,2],pch=19,cex=1)
# points(xyFromCell(S.tilde,as.numeric(names(n))),pch=19,cex=n/max(n)+0.25,col=2)
# y[10]

##################################################################################
### Simulate observed telemetry locations
##################################################################################

# Set parameters for observation mode
# lc <- s$lc  # location classes from PV95KOD09
lc <- sample(1:3,Ts,replace=TRUE)  # location classes from PV95KOD09
s <- get.s(mu.t,lc,sigma,a,rho,100)  # observed locations

# Plot true and observed locations
plot(s,pch="")
plot(S,add=TRUE)
points(mu.t,pch=19,cex=0.25,col=y+1) # All true locations
points(s,col=y+1,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu.t[,1],mu.t[,2],col="grey50") # Connections between s and mu
points(mu.t[y==1,],pch=19,col=rgb(1,1,1,0.6),cex=0.5) # Haul out locations

##################################################################################
### Write/read workspace and rasters
##################################################################################

# save.image('R/sim/sim.out.Rdata')
# load('R/sim/sim.out.Rdata')


##################################################################################
### Fit model
##################################################################################

I <- order(n,decreasing=TRUE)
start <- list(mu=as.numeric(names(n[I])),pi=n[I]/Ts,
	theta=theta,sigma.mu=sigma.mu,sigma=sigma,a=a,rho=rho,
	alpha=alpha,sigma.alpha=1,beta=beta) 

theta.priors <- DPelicit(Ts,mean=2,std=2,method="JGL")$inp  
hist(rgamma(1000,theta.priors[1],theta.priors[2]),breaks=100)
	# Gamma(r,q) prior for theta

priors <- list(mu.sigma=sigma.mu,sigma.sigma=0.25,f.hat=f.hat,
	sigma.beta=2,u.sigma=200000,r.alpha=2,q.alpha=1,
	r.theta=theta.priors[1],q.theta=theta.priors[2])
	# r.theta=1,q.theta=1)

tune <- list(mu=1290,sigma.mu=590,
	sigma=c(487,862,2629),rho=c(0.44,0.06,0.08),a=c(0.19,0.06,0.06))

source("~/Documents/research/harbor_seals/projects/central_places/R/central.places.mcmc.R")

out1 <- central.places.mcmc(s,lc,y,X.tmp,W,J=50,S.tilde.mat,
	priors=priors,tune=tune,start=start,adapt=TRUE,n.mcmc=5000)

out2 <- central.places.mcmc(s,lc,y,X.tmp,W,J=50,S.tilde.mat,
	priors=priors,tune=out1$tune,start=out1$start,adapt=FALSE,n.mcmc=20000)

save.image('R/sim/sim.out.Rdata')


##################################################################################
### Inspect model output
##################################################################################

mod <- out1 
mod <- out2

idx <- 1:1000
idx <- 1:5000
idx <- 1:10000
idx <- 1:20000

# Inference on mu(t): haul-out site locations
n.tmp <- table(mod$mu[,idx])  # tabulate posterior for mu.0
mu.t.post <- S.tilde
values(mu.t.post) <- NA
mu.t.post[as.numeric(names(n.tmp))] <- (n.tmp/max(n.tmp))^(1/2)
plot(mu.t.post,xlim=range(S.tilde.xy[,1]),ylim=range(S.tilde.xy[,2]))  # posterior of mu
points(xyFromCell(S.tilde,as.numeric(names(n))),pch=1,cex=n/max(n)+0.25,col=2)
points(s,pch=19,cex=0.2,col=3)

# Inference on mu_j: 'potential' haul-out site locations
n.tmp <- table(unlist(apply(mod$mu[,idx],2,unique)))
mu.j.post <- S.tilde
values(mu.j.post) <- NA
mu.j.post[as.numeric(names(n.tmp))] <- (n.tmp/max(n.tmp))^(1/1)
plot(mu.j.post,xlim=range(S.tilde.xy[,1]),ylim=range(S.tilde.xy[,2]))  # posterior of mu
points(xyFromCell(S.tilde,as.numeric(names(n))),pch=1,cex=n/max(n)+0.25,col=2)
points(s,pch=19,cex=0.2,col=3)

# Examine inference of particular telemetry observation
pt.idx <- 100
points(xyFromCell(S.tilde,mod$mu[pt.idx,idx]),pch=19,cex=0.5,col=rgb(0,0,0,0.025))
points(s.tmp[pt.idx,1],s.tmp[pt.idx,2],pch=19,col=5)

# Inference on beta: fixed effects of temporal haul-out process
matplot(mod$beta[idx,],type="l");abline(h=beta,col=1:5)
beta.hat <- apply(mod$beta[idx,],2,mean)
beta.quant <- t(apply(mod$beta[idx,],2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19,col=rgb(1,0,0,0.25))

# Inference on alpha: random effects of temporal haul-out process
alpha.hat <- apply(mod$alpha[idx,],2,mean)
alpha.quant <- t(apply(mod$alpha[idx,],2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
segments(1:qW,alpha.quant[,1],1:qW,alpha.quant[,2],col="lightgrey")
abline(h=0,col=2,lty=2)
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
points(alpha,pch=19,col=rgb(1,0,0,0.25))

# Inference on theta: DP concentration parameter
hist(mod$theta[idx],breaks=100)
plot(mod$theta[idx],type="l")
mean(mod$theta[idx])*log(Ts)

# Inference on m: modeled number of clusters
plot(mod$m,type="l")  # number of clusters  
barplot(table(mod$m)) 

# Inference on sigma_mu: homerange dispersion parameter
plot(mod$sigma.mu[idx],type="l")
hist(mod$sigma.mu[idx],breaks=100)
mean(mod$sigma.mu[idx])

# Inference on observation model parameters
matplot(mod$sigma,type="l");abline(h=sigma,col=1:3)
sigma.hat <- apply(mod$sigma[idx,],2,mean)
sigma.quant <- t(apply(mod$sigma[idx,],2,quantile,c(0.025,0.975)))
plot(sigma.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(sigma.quant)))
segments(1:3,sigma.quant[,1],1:3,sigma.quant[,2],col="lightgrey")
points(sigma.hat,pch=19,col=rgb(0,0,0,0.25))
points(sigma,pch=19,col=rgb(1,0,0,0.25))

matplot(mod$a,type="l");abline(h=a,col=1:3)
a.hat <- apply(mod$a[idx,],2,mean)
a.quant <- t(apply(mod$a[idx,],2,quantile,c(0.025,0.975)))
plot(a.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:3,a.quant[,1],1:3,a.quant[,2],col="lightgrey")
points(a.hat,pch=19,col=rgb(0,0,0,0.25))
points(a,pch=19,col=rgb(1,0,0,0.25))

matplot(mod$rho,type="l");abline(h=rho,col=1:3)
rho.hat <- apply(mod$rho[idx,],2,mean)
rho.quant <- t(apply(mod$rho[idx,],2,quantile,c(0.025,0.975)))
plot(rho.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:3,rho.quant[,1],1:3,rho.quant[,2],col="lightgrey")
points(rho.hat,pch=19,col=rgb(0,0,0,0.25))
points(rho,pch=19,col=rgb(1,0,0,0.25))

# Inference for v(t): auxilliary variable for haul-out process of telemetry locations
v.inv <- matrix(pnorm(mod$v[1:Ts,idx]),,Ts,byrow=TRUE)
v.inv.mean <- apply(v.inv,2,mean)
v.inv.quant <- t(apply(v.inv,2,quantile,c(0.025,0.975)))
plot(v.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:Ts,v.inv.quant[,1],1:Ts,v.inv.quant[,2],col=rgb(0,0,0,0.15))
abline(h=0.5,col=2,lty=2)
points(y,col=3,pch=19,cex=0.5)





# Assess autocorrelation in the v's
pred <- X.tmp%*%beta.hat + W%*%alpha.hat
hist(pred)
hist(v.inv.mean)
plot(abs(v.inv.mean-pred),type="l")

idx <- which(s$deploy_idx==deploy.idx)
tm <- strptime(X$dtime[idx],format="%Y-%m-%d %H:%M:%S",tz="UTC")

d.tm <- difftime(tm[-1],tm[-Ts],units="hours")
resid <- v.inv.mean-pred

resid <- v.inv[109,]-pred

plot(resid)
d.resid <- abs(resid[-1]-resid[-Ts])

plot(d.tm,d.resid)
abline(v=1)

d <- 1/as.matrix(dist(as.numeric(tm)))
diag(d) <- 0
d[is.infinite(d)] <- 0

library(spdep)
lw <- mat2listw(d) 
lwW <- nb2listw(lw$neighbours, glist=lw$weights, style="W") 
moran.test(c(v.inv.mean-pred), lwW, alternative="two.sided") 

d <- as.matrix(dist(as.numeric(tm)))
d.resid <- as.matrix(dist(c(v.inv.mean-pred)))
plot(c(d),c(d.resid),pch=19,cex=0.1,col=rgb(0,0,0,0.1))

