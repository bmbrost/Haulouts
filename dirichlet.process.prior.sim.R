###
### Simulation 1-dimensional Dirichlet process prior (Gelman et al. 2014, section 23.2)
###

a0 <- 0.1  # concentration parameter
K <- 100  # total number of clusters 
P0 <- seq(1,K,1)  # base probability measure; discrete uniform distribution 

# Generate data according to Chinese restaurant process
n <- 100
k <- numeric(n)

k[1] <- sample(P0,1)
for (i in 2:n) {
  tab <- table(k)
  idx <- as.numeric(names(tab))>0
  occ <- as.numeric(names(tab)[idx])
  p <- rep(a0/(a0+n-1),100)
  p[occ] <- tab[idx]/(a0+n-1)
  k[i] <- sample(P0,1,prob=p)
}

hist(k,breaks=0:n)


###
### Fit model
###

source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
   start=list(a0=a0),n.mcmc=1000)
tab <- hist(k,breaks=0:n,plot=FALSE)
boxplot(c(out1$P)*n~rep(1:n,out1$n.mcmc),pch=19,cex=0.25)
points(1:n,tab$counts,col="red")
points(1:n,apply(out1$P*n,1,mean),col="blue")

