###
### Simulation 1-dimensional Dirichlet process mixture
###

# Using blocked Gibbs sampler per Gelman et al. 2014, section 23.2
a0 <- 1  # concentration parameter
P0 <- seq(0,100,1)  # base probability measure; discrete uniform distribution 
n <- 100  # number of observations

# Generate data according to stick-breaking process

# Discrete uniform base distribution
N <- 20  # maximum number of clusters for truncation approximation to DPM
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))
theta <- sample(P0,N,replace=FALSE)  # clusters randomly drawn from P0
k <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments of observations
hist(k,breaks=P0)
tab.k <- table(k)

# Continuous uniform base distribution
N <- 20  # maximum number of clusters for truncation approximation to DPM
v <- c(rbeta(N-1,1,a0),1)
pie <- v*c(1,cumprod((1-v[-N])))

theta <- runif(N,min(P0),max(P0))  # clusters randomly drawn from P0
k <- sample(theta,n,replace=TRUE,prob=pie)  # cluster assignments of observations
hist(k,breaks=1000)
# hist(k,breaks=P0)  # for comparison to discrete base distribution
k.tab <- table(k)



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





# # Generate data according to Chinese restaurant process (Gelman et al. 2014, section 23.2)
# a0 <- 0.1  # concentration parameter
# K <- 100  # total number of clusters 
# P0 <- seq(1,K,1)  # base probability measure; discrete uniform distribution 
# N <- 20  # maximum number of clusters
# 
# n <- 100
# k <- numeric(n)
# 
# k[1] <- sample(P0,1)
# for (i in 2:n) {
#   tab <- table(k)
#   idx <- as.numeric(names(tab))>0
#   occ <- as.numeric(names(tab)[idx])
#   p <- rep(a0/(a0+n-1),100)
#   p[occ] <- tab[idx]/(a0+n-1)
#   k[i] <- sample(P0,1,prob=p)
# }
# 
# hist(k,breaks=0:n)
# 
# 
# ###
# ### Fit model
# ###
# 
# source("/Users/brost/Documents/git/Haulouts/dirichlet.process.prior.mcmc.R")
# out1 <- dirichlet.process.prior.mcmc(k,P0,priors=list(a=0.01,b=0.1),tune=list(a0=0.1),
#    start=list(a0=a0),n.mcmc=1000)
# tab <- hist(k,breaks=0:n,plot=FALSE)
# boxplot(c(out1$P)*n~rep(1:n,out1$n.mcmc),pch=19,cex=0.25)
# points(1:n,tab$counts,col="red")
# points(1:n,apply(out1$P*n,1,mean),col="blue")
# 
