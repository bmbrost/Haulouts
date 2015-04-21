###
### Simulate 1-dimensional haul out data 
###

# Define 1-dimensional space
L <- c(-10,-0.5) # Land domain
H <- c(-0.5,0.5) # Support of haul out process
S <- c(-0.5,10) # Support of movement process (marine and haul-out environments)

plot(0,0,xlim=c(-10,10),ylim=c(0,1),pch="",yaxt="n",xaxt="n",xlab="",ylab="")
polygon(x=c(L,rev(L),L[1]),y=c(0,0,1,1,1),col="gray45")
polygon(x=c(S,rev(S),S[1]),y=c(0,0,1,1,1),col="gray85")
polygon(x=c(H,rev(H),H[1]),y=c(0,0,1,1,1),angle=45,density=5)

n <- 100 # Number of locations
p <- 0.5 # Probability of being on the haul-out
z <- rbinom(n,1,p)
mu <- cbind(0,runif(n,0,1))
mu[z==1,1] <- runif(sum(z),H[1],H[2])
mu[z==0,1] <- runif(n-sum(z),S[1],S[2])

points(mu)
points(mu[z==1,],pch=19,col=rgb(1,1,1,0.6)) # Haul out locations


###
### Fit model
###

sum(z)

source("haulout.mcmc.R")
out1 <- haulout.mcmc(mu[,1],H,S,priors=list(alpha=1,beta=1),n.mcmc=10000)

plot(out1$p,type="l")
mean(out1$p)

warnings()


S
p/(p+p*dist(S))
dist(S)
1/11.5
?rbeta



