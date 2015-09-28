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
library(rgeos)
library(MCMCpack)  # for rdirichlet(...)
library(splines)
library(DPpackage)

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
### Simulate resource surface for haul-out RSF
##################################################################################

S.tilde.idx <- which(values(S.tilde)>0)
S.tilde.xy <- xyFromCell(S.tilde,S.tilde.idx)

U <- cbind(1,sqrt((S.tilde.xy[,1]-70000)^2+(S.tilde.xy[,2]-832000)^2))
qU <- ncol(U)

U.mean <- c(0,apply(U,2,mean)[-1])
U.sd <- c(1,apply(U,2,sd)[-1])
U.scale <- t(apply(U,1,function(x) (x-U.mean)/U.sd))

gamma <- c(1,-9000)  # coefficients on unstandarized U

# Examine RSF surface
selection <- exp(U.scale%*%(gamma/U.sd))  # RSF
S.tilde.rsf <- S.tilde
S.tilde.rsf[S.tilde.idx] <- selection
plot(S.tilde.rsf)
plot(U[,-1],selection,type="p")


##################################################################################
### Simulate haul-out sites and assignments using a stick-breaking process
### See Ishwaran and James (2001), Gelman et al. (2014), Section 23.2
##################################################################################

T <- 500  # number of locations to simulate
n <- 500  # number of wet/dry observations to simulate
theta <- 0.5  # Dirichlet process mixture concentration parameter
H <- 50  # maximum number of clusters for truncation approximation

E.m <- theta*(digamma(theta+T)-digamma(theta))  # expected number of clusters
theta.priors <- DPelicit(T,mean=E.m,std=3,method="JGL")$inp

# Simulate haul-out sits of telemetry locations s
S.tilde.idx <- which(values(S.tilde)>0)
mu.0 <- sample(S.tilde.idx,H,prob=selection)  # clusters randomly drawn from S.tilde
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

# Set parameters for observation model
sigma <- c(2291,2727,13252)  # standard deviation of error along x-axis 
# nu <- c(16.74,1.60,1.00)  # degrees of freedom for mixture t distribution
a <- c(0.70,0.50,0.75)  # modifies sd for error along y-axis
rho <- c(0.85,0.16,0.30)  # rotation in mixture components
lc <- sample(1:3,T,replace=TRUE)  # location classes
s <- get.s(mu,lc,sigma,a,rho,nu)  # observed locations


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
start <- list(theta=theta,ht=ht,z=z,pie=pie,beta=beta,gamma=gamma/U.sd,
	sigma=sigma,sigma.mu=sigma.mu,sigma.alpha=2) 
priors <- list(sigma.mu.l=0,sigma.mu.u=5000,sigma.beta=10,sigma.gamma=10,
	H=H,r.theta=theta.priors[1],q.theta=theta.priors[2],r.sigma.alpha=2,q.sigma.alpha=1,
	sigma=sigma,a=a,rho=rho,lc=lc)  # observation model parameters; empirical Bayes
tune <- list(mu.0=2000,sigma.mu=1250,gamma=2.0)
source("/Users/brost/Documents/git/haulouts/haulouts.3.mcmc.R")
out1 <- haulouts.3.mcmc(s,y,X=rbind(X[s.idx,],X[-s.idx,]),W=rbind(W[s.idx,],W[-s.idx,]),
	U.scale,S.tilde,priors=priors,tune=tune,start=start,n.mcmc=10000)


##################################################################################
### Inspect model output
##################################################################################

mod <- out1 
idx <- 1:1000
idx <- 1:2000
idx <- 1:5000
idx <- 1:10000

# Inference on haul-out site locations (mu.0)
tab.tmp <- table(mod$ht[,idx])
S.post <- S.tilde-1
S.post[as.numeric(names(tab.tmp))] <- (tab.tmp/max(tab.tmp))^(1/1)
plot(S.post)
points(xyFromCell(S.tilde,as.numeric(names(tab))),pch=1,cex=tab/max(tab)+0.25,col=2)
points(s,pch=19,cex=0.2,col=3)

pt.idx <- 1
points(xyFromCell(S.tilde,mod$ht[pt.idx,idx]),pch=19,cex=0.5,col=rgb(0,0,0,0.025))
points(mu[pt.idx,1],mu[pt.idx,2],pch=19,cex=0.75,col=5 )
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

# Haul-out resource selection covariates
matplot(mod$gamma[idx,-1],type="l"); abline(h=gamma[-1],col=1:(qU-1),lty=2) 
gamma.hat <- apply(mod$gamma[idx,],2,mean)
gamma.quant <- t(apply(mod$gamma[idx,],2,quantile,c(0.025,0.975)))
plot(gamma.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(gamma.quant)))
abline(h=0,col=2,lty=2)
segments(1:qU,gamma.quant[,1],1:qU,gamma.quant[,2],col="lightgrey")
points(gamma.hat,pch=19,col=rgb(0,0,0,0.25))
points(gamma/U.sd,pch=19)

# Standard deviation of random effects (sigma.alpha)
hist(mod$sigma.alpha[idx],breaks=100)

# Concentration parameter
hist(mod$theta[idx],breaks=100);abline(v=theta,col=2,lty=2) 
mean(mod$theta[idx])*log(T)

# Modeled number of clusters
plot(mod$m,type="l");abline(h=m,col=2,lty=2)  # true number of clusters  
barplot(table(mod$m)) 

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
