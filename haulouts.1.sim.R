##################################################################################
### Simulate 2-dimensional haul out data using a Dirichlet process mixture 
### Developed for application with real harbor seal data using rasters, 
### mixture t distribution, etc.
##################################################################################

rm(list=ls())

# library(cluster)
# library(msm)  # for simulating mu
# library(pscl)
library(sp)
library(rgdal)
library(maptools)
library(raster)
# library(gstat)
# library(rgeos)
library(MCMCpack)  # for rdirichlet(...)
library(splines)


##################################################################################
### Prepare spatial data for analysis
##################################################################################

###
### Define S.tilde (support of haul-out sites, mu_0)
###

proj.aea <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 
	+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
ak.simple <- readOGR(dsn=
	"/Users/brost/Documents/research/harbor_seals/geodata/base/physical_features",
	layer="ak_simple_poly_R")
ak <- readOGR(dsn=
	"/Users/brost/Documents/research/harbor_seals/geodata/base/physical_features", 
	layer="ak_63360_poly_R")

# Define S according to some bounding box
S.clip <- matrix(c(55000,110000,110000,55000,55000,805000,805000,840000,840000,805000),5,2)
S.clip <- SpatialPolygons(list(Polygons(list(Polygon(S.clip)),"1")),proj4string=CRS(proj.aea))
S.poly <- gIntersection(ak.simple,S.clip)

# Create raster of S
S.res <- 250
S <- raster(S.poly,resolution=S.res)
S <- rasterize(S.poly,S,field=1)
S <- reclassify(S,matrix(c(1,NA,NA,1),2,2,byrow=TRUE))
plot(S)

# Identify cells along shoreline to define S.tilde
S.tilde <- raster(S.poly,resolution=S.res)
S.tilde <- rasterize(S.poly,S.tilde,field=1)
S.tilde <- boundaries(S.tilde,asNA=TRUE)
plot(S.tilde)
plot(S,colNA=NA,add=TRUE,col="grey85")
plot(S.poly,add=TRUE)


##################################################################################
### Simulate haul-out sites and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
##################################################################################

T <- 500  # number of locations to simulate
n <- 500  # number of wet/dry observations to simulate
theta <- 2.0  # Dirichlet process mixture concentration parameter
H <- 50  # maximum number of clusters for truncation approximation

# Simulate haul-out sits of telemetry locations s
S.tilde.idx <- which(values(S.tilde)>0)
mu.0 <- sample(S.tilde.idx,H)  # clusters randomly drawn from S.tilde
plot(S.tilde)
points(xyFromCell(S.tilde,mu.0),col=2,pch=19,cex=0.25)  # plot cluster locations
eta <- c(rbeta(H-1,1,theta),1)  # stick-breaking weights
pie <- eta*c(1,cumprod((1-eta[-H])))  # probability mass
plot(pie,type="b")
ht <- sample(mu.0,T,replace=TRUE,prob=pie)  # latent cluster assignments
m <- length(unique(ht))
tab <- table(ht)
plot(S.tilde)
points(xyFromCell(S.tilde,as.numeric(names(tab))),pch=19,cex=tab/max(tab)+0.25,col=2)


##################################################################################
### Simulate wet/dry status for T telemetry locations (s) and T SEA records (y)
##################################################################################

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
y <- rbinom(n,1,p[-s.idx])  # Haul-out indicator variable for y: 1=hauled-out, 0=at-sea
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


##################################################################################
### Simulate true and observed telemetry locations
##################################################################################

# Simulate true locations
sigma.mu <- 2000  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)

# Sample true locations for hauled-out status; i.e., true location=haul-out site
mu[z==1,] <- xyFromCell(S.tilde,ht[z==1])

# Sample true locations for at-sea locations
idx <- which(values(S)==1)
S.xy <- xyFromCell(S,idx)
mu.0.xy <- xyFromCell(S.tilde,ht[z==0])
mu[z==0,1] <- apply(mu.0.xy,1,function(x) sample(idx,1,prob=
	dnorm(x[1],S.xy[,1],sigma.mu)*dnorm(x[2],S.xy[,2],sigma.mu))) 
mu[z==0,] <- xyFromCell(S,mu[z==0,1])


# Check true location generation
mu.0.xy <- xyFromCell(S.tilde,ht)
S.test <- S
S.test[idx] <- dnorm(mu.0.xy[14,1],S.xy[,1],sigma.mu)*dnorm(mu.0.xy[14,2],S.xy[,2],sigma.mu)
plot(S.tilde)
plot(S.test,add=TRUE)
points(mu.0.xy[14,1],mu.0.xy[14,2])
points(mu,pch=19,cex=0.25)
points(xyFromCell(S.tilde,as.numeric(names(tab))),pch=19,cex=tab/max(tab)+0.25,col=2)


# Simulate observed locations
sigma <- 5000 # Observation error
s <- mu
s <- s+rnorm(T*2,0,sigma) # Add error to true locations

# Plot true and observed locations
plot(S.tilde,col=rgb(1,0,0))  # plot support
plot(S,add=TRUE)
points(mu,pch=19,cex=0.25,col=z+1) # All true locations
points(s,col=z+1,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6),cex=0.5) # Haul out locations


##################################################################################
### Fit models
##################################################################################

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulouts.1.mcmc.R")
start <- list(theta=theta,ht=ht,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie,beta=beta)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r=2,q=0.1,sigma.l=0,sigma.u=10000,sigma.mu.l=0,sigma.mu.u=5000,
	sigma.beta=10)
tune <- list(mu.0=3500,sigma=750,sigma.mu=1500)
# hist(rgamma(1000,2,0.1))
# hist(rgamma(1000,50,10))
out1 <- haulouts.1.mcmc(s,y,X[s.idx,],X[-s.idx,],W[s.idx,],W[-s.idx,],
	S.tilde,sigma.alpha=2,priors=priors,tune=tune,start=start,n.mcmc=2000)
hist(rbeta(1000,1,4))

##################################################################################
### Inspect model output
##################################################################################

mod <- out1 
idx <- 1:1000
idx <- 1:2000
idx <- 1:10000

# Inference on haul-out site locations (mu.0)
tab.tmp <- table(mod$ht[,idx])
S.post <- S.tilde-1
S.post[as.numeric(names(tab.tmp))] <- (tab.tmp/max(tab.tmp))^(1/2)
plot(S.post)
points(xyFromCell(S.tilde,as.numeric(names(tab))),pch=1,cex=tab/max(tab)+0.25,col=2)
points(s,pch=19,cex=0.2,col=3)

pt.idx <- 100
points(xyFromCell(S.tilde,mod$ht[pt.idx,idx]),pch=19,cex=0.5,col=rgb(0,0,0,0.025))
points(mu[pt.idx,1],mu[pt.idx,2],pch=19,cex=0.75,col=5)
points(s[pt.idx,1],s[pt.idx,2],pch=19,col=2)
table(mod$z[pt.idx,idx])
z[pt.idx]

# Haul-out probability covariates
matplot(mod$beta[idx,],type="l"); abline(h=beta,col=1:3,lty=2) 
beta.hat <- apply(mod$beta[idx,],2,mean)
beta.quant <- t(apply(mod$beta[idx,],2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))
points(beta,pch=19)

# Inference on alpha
alpha.hat <- apply(mod$alpha[idx,],2,mean)
alpha.quant <- t(apply(mod$alpha[idx,],2,quantile,c(0.025,0.975)))
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
hist(mod$theta[idx],breaks=100);abline(v=theta,col=2,lty=2) 
mean(mod$theta[idx])*log(T)

# Modeled number of clusters
plot(mod$m,type="l");abline(h=m,col=2,lty=2)  # true number of clusters  
barplot(table(mod$m)) 

# Observation error
hist(mod$sigma[idx],breaks=100);abline(v=sigma,col=2,lty=2)
mean(mod$sigma[idx])

# Dispersion about haul-out for at-sea locations
hist(mod$sigma.mu[idx],breaks=100);abline(v=sigma.mu,col=2,lty=2)
mean(mod$sigma.mu[idx])

# Inference on z: latent haul-out indicator variable for telemetry locations
# Note: estimation of z (z.hat) is based on covariates and location of telemetry 
# observations, whereas u is based on covariates alone.
z.hat <- apply(mod$z[,idx],1,sum)/(length(idx))
boxplot(z.hat~z,ylim=c(0,1))

v <- apply(mod$beta[idx,],1,function(x) X[s.idx,]%*%x)+
	apply(mod$alpha[idx,],1,function(x) W[s.idx,]%*%x)
v.inv <- matrix(pnorm(v),,T,byrow=TRUE)
v.inv.mean <- apply(v.inv,2,mean)
v.inv.quant <- t(apply(v.inv,2,quantile,c(0.025,0.975)))
plot(v.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,v.inv.quant[,1],1:T,v.inv.quant[,2],col=rgb(0,0,0,0.15))
abline(h=0.5,col=2,lty=2)
points(z,col=3,pch=19,cex=0.5)
points(z.hat,pch=19,cex=0.25)

# Inference for y
v.tilde.inv <- matrix(pnorm(mod$v[-(1:T),idx]),,n,byrow=TRUE)
v.tilde.inv.mean <- apply(v.tilde.inv,2,mean)
v.tilde.inv.quant <- t(apply(v.tilde.inv,2,quantile,c(0.025,0.975)))
plot(v.tilde.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:n,v.tilde.inv.quant[,1],1:n,v.tilde.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y,col=4,pch=19,cex=0.5)
