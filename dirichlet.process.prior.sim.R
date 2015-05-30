###
### Simulation 1-dimensional Dirichlet process prior (Gelman et al. 2014, section 23.2)
###

a0 <- 0.1  # concentration parameter
P0 <- seq(0,6,0.0001)  # base probability measure
n <- 100  # number of observations

# Generate data according to Chinese restaurant process

# Discrete uniform base distribution
k <- numeric(n)
k[1] <- sample(P0,1)
for (i in 2:n) {  
  tab <- table(k)
  idx <- as.numeric(names(tab))>0
  occ <- as.numeric(names(tab)[idx])+1
  p <- rep(a0/(a0+n-1),length(P0))
  p[occ] <- tab[idx]/(a0+n-1)
  k[i] <- sample(P0,1,prob=p)
}
hist(k,breaks=1000)
k.tab <- table(k)

# Continuous uniform base distribution
k <- numeric(n)
k[1] <- runif(1,min(P0),max(P0))
for (i in 2:n) {
  tab <- table(k[k>0])
  p <- tab/(a0+n-1)
  p <- c(p,1-sum(p))
  k.idx <- sample(1:length(p),1,prob=p)
  k[i] <- ifelse(k.idx==length(p),runif(1,min(P0),max(P0)),as.numeric(names(tab)[k.idx]))
}
hist(k,breaks=1000)
k.tab <- table(k)


###
### Fit model
###

source("/Users/bmb/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
   start=list(a0=a0),n.mcmc=1000)
# boxplot(c(out1$P)*n~rep(as.numeric(names(k.tab)),out1$n.mcmc),pch=19,cex=0.25)
plot(rep(as.numeric(names(k.tab)),out1$n.mcmc),c(out1$P)*n,pch=19,cex=0.25,col=rgb(0,0,0,0.05))
points(as.numeric(names(k.tab)),k.tab,col="red")
points(as.numeric(names(k.tab)),apply(out1$P*n,1,mean),col="blue")

