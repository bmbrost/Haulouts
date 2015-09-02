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


#############################################################################################
### Prepare spatial data for analysis
#############################################################################################

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

# Define S according to centroid
# center <- matrix(c(61156,830265),1) #PV95KOD12
# center <- SpatialPoints(center,proj4string=CRS(proj.aea))

# Create polygon of S
# S.buffer <- gBuffer(center,width=20000)
# S.poly <- gIntersection(ak,S.buffer,byid=TRUE)
# plot(S.poly,col="grey")

# Create raster of S
S.res <- 250
S <- raster(S.poly,resolution=S.res)
S <- rasterize(S.poly,S,field=1)
S <- reclassify(S,matrix(c(1,NA,NA,1),2,2,byrow=TRUE))
plot(S)

# Remove isolated water cells enforcing 4-neighbor rule
# Note: make sure haulout is inside S
# tran <- transition(S,function(x) 1,directions=4)
# idx <- which(values(S)==1)
# cd <- costDistance(tran, haulout, xyFromCell(S,idx))*xres(S)
# values(S)[idx[which(cd==Inf)]] <- NA

# Limit extent of S to some distance d from haul-out location
# d <- 60000
# values(S)[idx[which(cd>d)]] <- NA
# plot(S)

# Identify cells along shoreline to define S.tilde
S.tilde <- raster(S.poly,resolution=S.res)
S.tilde <- rasterize(S.poly,S.tilde,field=1)
S.tilde <- boundaries(S.tilde,asNA=TRUE)
# idx <- which(values(S.tilde)==1)
# S.tilde[idx] <- 1:length(idx)  # values of haul-out cells labeled in sequence
plot(S.tilde)
plot(S,colNA=NA,add=TRUE,col="grey85")
plot(S.poly,add=TRUE)

# Convert to points
# idx <- which(values(S.tilde)==1)
# S.tilde.pts <- xyFromCell(S.tilde,idx)
# S.tilde.pts <- data.frame(S.tilde.pts,survy_d="shoreline",depth=0,qlty_cd=NA,active=NA)
# S.tilde.pts <- SpatialPointsDataFrame(S.tilde.pts[,1:2], S.tilde.pts[,3:6],
	# proj4string=CRS(proj.aea))
# points(S.tilde.pts,pch=19,cex=0.1)


#############################################################################################
### Simulate haul-out sites and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
#############################################################################################

T <- 500  # number of locations to simulate
n <- 500  # number of wet/dry observations to simulate
theta <- 1.0  # Dirichlet process mixture concentration parameter
H <- 25  # maximum number of clusters for truncation approximation

# Simulate haul-out sits of telemetry locations s
idx <- which(values(S.tilde)>0)
# points(xyFromCell(S.tilde,idx))
mu.0 <- sample(idx,H)  # clusters randomly drawn from S.tilde
mu.0 <- xyFromCell(S.tilde,mu.0)  # clusters locations
points(mu.0,col=1,pch=19,cex=0.25)  # plot cluster locations
eta <- c(rbeta(H-1,1,theta),1)  # stick-breaking weights
pie <- eta*c(1,cumprod((1-eta[-H])))  # probability mass
h.match <- sample(1:H,T,replace=TRUE,prob=pie)  # latent cluster assignments
h <- mu.0[h.match,]  # latent clusters
points(mu.0[unique(h.match),],pch=19,cex=table(h.match)/max(table(h.match))+0.25,col=2)


###
### Simulate wet/dry status for T telemetry locations (s) and T SEA records (y)
###

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


###
### Simulate true and observed locations
###

# Simulate true locations
sigma.mu <- 2000  # dispersion about haul-out for at-sea locations
mu <- matrix(0,T,2)
mu[z==1,] <- h[z==1,]  # true location is haul-out site for hauled-out locations

# Sample true locations for at-sea locations
idx <- which(values(S)==1)
S.xy <- xyFromCell(S,idx)
mu[z==0,1] <- apply(h[z==0,],1,function(x) sample(idx,1,prob=
	dnorm(x[1],S.xy[,1],sigma.mu)*dnorm(x[2],S.xy[,2],sigma.mu))) 
mu[z==0,] <- xyFromCell(S,mu[z==0,1])

# Check true location generation
S.test <- S
S.test[idx] <- dnorm(mu.0[3,1],S.xy[,1],sigma.mu)*dnorm(mu.0[3,2],S.xy[,2],sigma.mu)
plot(S.tilde)
plot(S.test,add=TRUE)
points(mu,pch=19,cex=0.25)
points(mu.0[unique(h.match),],pch=19,cex=table(h.match)/max(table(h.match))+0.25,col=2)

# Simulate observed locations
sigma <- 5000 # Observation error
s <- mu
s <- s+rnorm(T*2,0,sigma) # Add error to true locations


# Plot support
plot(S.tilde)
plot(S,add=TRUE)

# Plot true and observed locations
points(mu,pch=19,cex=0.25,col=z+1) # All true locations
points(s,col=z+1,pch=3,cex=0.5) # Observed locations
segments(s[,1],s[,2],mu[,1],mu[,2],col="grey50") # Connections between s and mu
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6),cex=0.5) # Haul out locations


###
### Fit models
###

# Fit model using blocked Gibbs sampler 
source("/Users/brost/Documents/git/haulouts/haulouts.1.mcmc.R")
start <- list(theta=theta,h=h,z=z,p=p,#h=fitted(kmeans(s,rpois(1,10))),
  sigma=sigma,sigma.mu=sigma.mu,pie=pie,beta=beta)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r=4,q=2,sigma.l=0,sigma.u=10000,sigma.mu.l=0,sigma.mu.u=5000,sigma.beta=1)
tune <- list(mu.0=1500,sigma=300,sigma.mu=100)
# hist(rgamma(1000,4,2))
# hist(rgamma(1000,5,2.5))
out1 <- haulouts.1.mcmc(s,y,X[s.idx,],X[-s.idx,],W[s.idx,],W[-s.idx,],
	S.tilde,sigma.alpha=2,priors=priors,tune=tune,start=start,n.mcmc=1000)

mod <- out1 
# mod <- out2
# mod <- out3
# mod <- out4
# dev.new()
# idx <- 1:100
idx <- 1:1000
idx <- 1:2000
idx <- 1:20000

# True clusters
plot(S.tilde)
plot(S,add=TRUE)
str(mod)

points(mod$h[,1,idx],mod$h[,2,idx],pch=19,cex=0.5,col=rgb(0,0,0,0.002))
points(h,pch=1,cex=1,col=rgb(1,0,0,1))
points(s,pch=19,cex=0.2,col=3)

pt.idx <- 137
points(mod$h[pt.idx,1,idx],mod$h[pt.idx,2,idx],pch=19,cex=0.2,col=rgb(1,1,0,0.25))
points(mu[pt.idx,1],mu[pt.idx,2],pch=19,cex=0.75)
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
hist(mod$a0[idx],breaks=100);abline(v=a0,col=2,lty=2) 
mean(mod$a0[idx])*log(T)

# Modeled number of clusters
# n.cls <- apply(mod$h[,,idx],3,function(x) nrow(unique(x)))
plot(mod$n.cls,type="l");abline(h=nrow(unique(h)),col=2,lty=2)  # true number of clusters  
barplot(table(mod$n.cls)) 

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
plot(s[,1],z.hat,ylim=c(0,1),col=z+1)

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
