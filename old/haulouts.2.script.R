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
library(mvtnorm)
library(rgeos)
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
proj.longlat <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
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
### Load telemetry data and subset
##################################################################################

#KODIAK SEALS
ko <- read.csv("~/Documents/research/harbor_seals/data/SDR data_CSU_v3/Argos/Kodiak 1993-1996-Argos.csv")
names(ko)[1] <- "deploy.idx"
head(ko)
ko$dtime <- strptime(ko$DateTime, format="%m/%d/%Y %H:%M",tz="GMT")  # in GMT/UTC
ko$lc <- factor(ko$LocationQuality, levels=c("3","2","1","0","A","B"), ordered=TRUE)
head(ko)
dup.idx <- duplicated(ko[,c(1,4)])
table(dup.idx)

# Subset seals
cbind(table(ko$deploy.idx,ko$lc),table(ko$deploy.idx))
apply(table(ko$deploy.idx,ko$Ptt),1,function(x) length(which(x==0)))
cbind(table(ko$Ptt,ko$lc),table(ko$deploy.idx))

# idx <- "PV95KOD09"
idx <- "PV95KOD12"
ko <- ko[which(ko$deploy.idx %in% idx),]
ko <- ko[-c(52,1306,1448),]  # exclude extreme points


# Define projections/Reproject
s <- ko
s.longlat <- SpatialPointsDataFrame(s[,9:8], s[,c(1,2,10:11)], proj4string=CRS(proj.longlat))
s <- spTransform(s.longlat, CRS=CRS(proj.aea)) 
s$deploy.idx <- factor(s$deploy.idx, levels=idx)
T <- nrow(s)
head(s)
range(s$dtime)
points(s@coords)

# Convert date-time to local AK time
# The Argos tag information document indicates that local time corresponding to 00:00 GMT
# for this inidividual is 13:00; therefore, subtract 11 hours to get AK time. It's not clear
# why this is the case because AKST = GMT-9 hours and AKDT = GMT-8 hours.

strptime("10/10/1995 00:00:00", format="%m/%d/%Y %H:%M:%S",tz="GMT")-(60*60*11)  # example
difftime(s$dtime-(11*60*60),s$dtime,units="hours")

s$dtime$hour <- s$dtime$hour-11  # convert to local time


##################################################################################
### Load SEA wet/dry status information
##################################################################################

# List of *.sea files for Kodiak seals
files <- list.files("~/Documents/research/harbor_seals/data/SDR data_CSU_v3/sea/raw/Kodiak_SEA",pattern="\\.SEA$",	full.names=TRUE)
length(files)

# Receptacle for records from all *.sea files
sea <- data.frame("Ptt"=NA,"Date"=NA,"Time"=NA,"Sea"=NA)[-1,]

# Number of lines to skip in each file
skip <- c(3,0,0,0,0,3,3,0,0,0,3,0,0,0)

n.sea <- 0 # Counter
for(i in 1:length(files)){
	sea.tmp <- read.csv(files[i],skip=skip[i],header=FALSE)
	if(ncol(sea.tmp)==5) sea.tmp <- sea.tmp[,-1] #Remove AnimalID column 
		# if present because it's kooky
	names(sea.tmp) <- names(sea)
	n.sea <- n.sea+nrow(sea.tmp)
	sea <-	rbind(sea,sea.tmp)
}

sea$Time <- formatC(as.numeric(sea$Time),format="d",width=6,flag="0")
sea$dtime <- paste(sea$Date,sea$Time)
sea$dtime <- strptime(sea$dtime, format="%m/%d/%Y %H%M%S",tz="GMT")  # in GMT/UTC
head(sea)

tapply(sea$dtime,sea$Ptt,range)

# Identify duplicate sea records
table(duplicated(sea))
idx <- !duplicated(sea)
dups <- sea[!idx,]

# Confirm duplicate records
sea[which(sea$dtime==dups[50,5]),]

# Eliminate duplicate records
sea <- sea[idx,]
table(duplicated(sea))

# Subset for individual seal
idx <- which(sea$Ptt==unique(s$Ptt))
y <- sea[idx,]
str(y)

# sea records for ptt 5046 span two time periods; take subset overlapping Argos locations
idx <- y$dtime>strptime("08/01/1995", format="%m/%d/%Y",tz="GMT")
range(y$dtime[idx])
range(y$dtime[!idx])
range(s$dtime)

y <- y[idx,]  # sea records corresponding to Argos locations for ptt 5046
range(y$dtime)
y$dry <- ifelse(y$Sea==1,0,1)
head(y);tail(y)

y$dtime$hour <- y$dtime$hour-11  # convert to local time

# write.csv(comb,"~/Documents/research/harbor_seals/data/SDR data_CSU_v3/sea/ko_sea.csv",row.names=FALSE)

##################################################################################
### Compare times at which Argos and Sea data were collected
##################################################################################

t.min <- sapply(1:nrow(s),function(x) which.min(abs(difftime(s$dtime[x],y$dtime))))
hist(as.numeric(difftime(s$dtime,y$dtime[t.min],units="hours")),breaks=1000)

t.tab <- sapply(1:nrow(s),function(x) sum(abs(difftime(s$dtime[x],y$dtime))<10))
hist(t.tab,breaks=20)

t.tab.idx <- sapply(1:nrow(s),function(x) which(abs(difftime(s$dtime[x],y$dtime))<10))
t.tab <- unlist(lapply(t.tab.idx,length))
hist(t.tab,breaks=20)

t.tab.sum <- lapply(t.tab.idx,function(x) sum(y$dry[x]))

hist(unlist(t.tab.sum)/t.tab,breaks=20)

plot(y$dtime,y$dry,ylim=c(0,1))
abline(v=as.numeric(s$dtime),col=2,lty=2)

idx.tmp <- 1000:1100
plot(y$dtime[idx.tmp],y$dry[idx.tmp],ylim=c(0,1))
abline(v=as.numeric(s$dtime),col=2,lty=2)



##################################################################################
### Define design matrix for haul-out use covariates
##################################################################################

X <- cbind(1,s$dtime$hour,s$dtime$yday)  # hour and day of year for Argos locations
X.tilde <- cbind(1,y$dtime$hour,y$dtime$yday)  # hour and day of year for SEA observations
# month(s$dtime)
qX <- ncol(X)

# Combine into single matrix
X <- rbind(X,X.tilde)  # Argos and SEA covariates
plot(X[,3])
# Transform hour-of-day to hours since midnight
X[,2] <- sapply(X[,2],function(x) min(24-x,x))

# Transform Julian day to number of days since 09 OCT 1995
dayone <- min(s$dtime,y$dtime)$yday  # Julian day corresponding to 09 OCT 1995
plot(cbind(X[,3],ifelse(X[,3]>=dayone,X[,3]-dayone,X[,3]+365-dayone)))
X[,3] <- ifelse(X[,3]>=dayone,X[,3]-dayone,X[,3]+365-dayone)

# Center and scale design matricies
X.mid <- apply(X[,-1],2,mean)
# X.scale <- apply(X[,-1],2, sd)
X.scale <- apply(X[,-1],2,function(x) max(x)-min(x))/6
X[,-1] <- t(apply(X[,-1],1,function(x) (x-X.mid)/X.scale))
summary(X)


##################################################################################
### Define B-spline basis expansion over time
##################################################################################

# Knots defined over time interval in which s and y occur
knots <- seq.POSIXt(min(s$dtime,y$dtime),max(s$dtime,y$dtime),by=60*60*2)  # 2-hour intervals

# knots <- seq(min(X[,3]),max(X[,3]),by=1/(X.scale[2]*12))  # knots at 2-hour intervals
knots <- knots[-c(1,length(knots))]  # remove knots at boudnaries

W <- bs(c(s$dtime,y$dtime),knots=knots,degree=3,,intercept=FALSE)  # cubic spline

matplot(W[1:nrow(s),],type="l",col=2,lty=1)
matplot(W[(nrow(s)+1):nrow(W),],type="l",col=2,lty=1)
# matplot(W[,],type="l",col=2,lty=1)


##################################################################################
### Get starting values
##################################################################################

library(dplyr)  # for dense_rank()

plot(S)
points(s)

cls <- locator()
cls <- cbind(cls$x,cls$y)

idx <- which(values(S.tilde)==1)
xy <- xyFromCell(S.tilde,idx)
mu.0 <- apply(cls,1,function(x) which.min(sqrt((x[1]-xy[,1])^2+(x[2]-xy[,2])^2)))
mu.0 <- idx[mu.0]
mu.0.xy <- xyFromCell(S.tilde,mu.0)

plot(S.tilde)
points(mu.0.xy,pch=19,cex=0.5)

ht <- apply(s@coords,1,function(x) which.min(sqrt((x[1]-mu.0.xy[,1])^2+(x[2]-mu.0.xy[,2])^2)))
ht <- mu.0[ht]

plot(S)
points(s,col=dense_rank(ht))

theta <- 1
z <- rbinom(T,1,0.5)
sigma.mu <- 5000

# Estimate normal distribution parameters from t-distribution parameters
sigma.t <- c(2259,1089,1365,2720,2702,13338)
nu <- c(16.52,2.29,1.72,1.59,1.23,1.01)
sigma <- sapply(1:6,function(x) sd(c(rmvt(10000000,sigma=sigma.t[x]*diag(2),df=nu[x]))))
a <- c(0.69,0.42,0.64,0.50,0.90,0.74)
rho <- c(0.85,0.73,0.40,0.16,0.21,0.30)
H <- 20
pie <- rdirichlet(1,rep(1/H,H))
beta <- c(-2.0,0.1,0.25)


##################################################################################
### Fit model
##################################################################################

theta.priors <- DPelicit(T,mean=3,std=3,method="JGL")$inp

# Fit model using blocked Gibbs sampler 
start <- list(theta=theta,ht=ht,z=z,pie=pie,beta=beta, 
  sigma.mu=sigma.mu,sigma.alpha=5)  # rdirichlet(1,rep(1/H,H))) 
priors <- list(H=H,r.theta=theta.priors[1],q.theta=theta.priors[2],
	r.sigma.alpha=2,q.sigma.alpha=1,sigma=sigma,a=a,rho=rho,lc=lc,
	sigma.l=0,sigma.u=10000,sigma.mu.l=0,sigma.mu.u=5000,sigma.beta=10)
tune <- list(mu.0=3500,sigma.mu=500)
# hist(rgamma(1000,2,4))
# hist(rgamma(1000,50,10))
source("/Users/brost/Documents/git/haulouts/haulouts.2.mcmc.R")
out1 <- haulouts.2.mcmc(s@coords,y$dry,X,W,S.tilde,
	priors=priors,tune=tune,start=start,n.mcmc=1000)

hist(out1$theta,breaks=100);abline(v=theta,col=2,lty=2) 
plot(out1$m,type="l");abline(h=m,col=2,lty=2)  # true number of clusters  
barplot(table(out1$m)) 
mean(out1$m)
sd(out1$m)
m

theta*(digamma(theta+T)-digamma(theta))  # expected number of clusters

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
points(s,pch=19,cex=0.1,col=3)
plot(S.post,add=TRUE)

pt.idx <- 59
points(xyFromCell(S.tilde,mod$ht[pt.idx,idx]),pch=19,cex=0.5,col=rgb(0,0,0,0.025))
points(s[pt.idx,1],s[pt.idx,2],pch=19,col=2)
table(mod$z[pt.idx,idx])

# Haul-out probability covariates
matplot(mod$beta[idx,],type="l")
beta.hat <- apply(mod$beta[idx,],2,mean)
beta.quant <- t(apply(mod$beta[idx,],2,quantile,c(0.025,0.975)))
plot(beta.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(beta.quant)))
abline(h=0,col=2,lty=2)
segments(1:qX,beta.quant[,1],1:qX,beta.quant[,2],col="lightgrey")
points(beta.hat,pch=19,col=rgb(0,0,0,0.25))

# Inference on alpha
alpha.hat <- apply(mod$alpha[idx,],2,mean)
alpha.quant <- t(apply(mod$alpha[idx,],2,quantile,c(0.025,0.975)))
plot(alpha.hat,pch=19,col=rgb(0,0,0,0.25),ylim=c(range(alpha.quant)))
abline(h=0,col=2,lty=2)
segments(1:ncol(W),alpha.quant[,1],1:ncol(W),alpha.quant[,2],col="lightgrey")
points(alpha.hat,pch=19,col=rgb(0,0,0,0.25))
# points(alpha,pch=19,col=3)

# Concentration parameter
hist(mod$theta[idx],breaks=100)
mean(mod$theta[idx])*log(T)

# Modeled number of clusters
plot(mod$m,type="l")
barplot(table(mod$m)) 

# Dispersion about haul-out for at-sea locations
hist(mod$sigma.mu[idx],breaks=100)
mean(mod$sigma.mu[idx])

# Inference on z: latent haul-out indicator variable for telemetry locations
# Note: estimation of z (z.hat) is based on covariates and location of telemetry 
# observations, whereas u is based on covariates alone.
z.hat <- apply(mod$z[,idx],1,sum)/(length(idx))
v <- apply(mod$beta[idx,],1,function(x) X[1:T,]%*%x)+
	apply(mod$alpha[idx,],1,function(x) W[1:T,]%*%x)
v.inv <- matrix(pnorm(v),,T,byrow=TRUE)
v.inv.mean <- apply(v.inv,2,mean)
v.inv.quant <- t(apply(v.inv,2,quantile,c(0.025,0.975)))
plot(v.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:T,v.inv.quant[,1],1:T,v.inv.quant[,2],col=rgb(0,0,0,0.15))
abline(h=0.5,col=2,lty=2)
points(z.hat,pch=19,cex=0.25)


t.min <- sapply(1:nrow(s),function(x) which.min(abs(difftime(s$dtime[x],y$dtime))))
hist(as.numeric(difftime(s$dtime,y$dtime[t.min],units="hours")))
data.frame(s$dtime,y$dtime[t.min])

which.max(difftime(s$dtime,y$dtime[t.min],units="hours"))

difftime(s$dtime[200],y$dtime[t.min[200]],units="hours")

idx.tmp <- 300:500
plot(y$dtime[idx.tmp],y$dry[idx.tmp],ylim=c(0,1))
abline(v=as.numeric(s$dtime[59]),col=2,lty=2)
abline(v=as.numeric(s$dtime),col=2,lty=2)
points(s$dtime,z.hat,pch=19,cex=0.5)

points(seq.POSIXt(min(s$dtime,y$dtime),max(s$dtime,y$dtime),by=60*60*2),alpha.hat[-1]/(max(abs(alpha.hat[-1]))*2)+0.5,lty=1,pch=19,col="gray50")

points(seq.POSIXt(min(s$dtime,y$dtime),max(s$dtime,y$dtime),by=60*60*2),alpha.hat[-1]/(max(abs(alpha.hat[-1]))*2)+0.5,lty=1,pch=19,col="gray50")

ncol(W)
W%*%alpha.hat

points(attr(W,"knots"),alpha.hat)
points(attr(W,"knots"),rep(0.5,ncol(W)-3))

ncol(mod$alpha)
length(alpha.hat)
length(seq.POSIXt(min(s$dtime,y$dtime),max(s$dtime,y$dtime),by=60*60*2)[1])
attr(W,"Boundary.knots")
attr(W,"knots")
ncol(W)

alpha.hat
points(as.numeric(s$dtime[20]),0.5,pch=19,cex=1,col=3)


T









# Inference for y
v.tilde.inv <- matrix(pnorm(mod$v[-(1:T),idx]),,nrow(y),byrow=TRUE)
v.tilde.inv.mean <- apply(v.tilde.inv,2,mean)
v.tilde.inv.quant <- t(apply(v.tilde.inv,2,quantile,c(0.025,0.975)))
plot(v.tilde.inv.mean,pch=19,col=rgb(0,0,0,0.25),ylim=c(0,1))
segments(1:nrow(y),v.tilde.inv.quant[,1],1:nrow(y),v.tilde.inv.quant[,2],col=rgb(0,0,0,0.15))
points(y$dry,col=4,pch=19,cex=0.5)
