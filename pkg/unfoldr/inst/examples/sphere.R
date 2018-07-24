library(unfoldr)

## beta distribution for radii
lam <- 10
theta <- list("size"=list("shape1"=2,"shape2"=4))

#debug(simPoissonSystem)
box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
S <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",pl=101)

# planar vertical section
sp <- planarSection(S,d=2.5)
## or equivalently
## sp0 <- 2.0*sphereIntersection(S,2.5,c(0,0,1),FALSE,10)
## sp[1:10]
## sp0[1:10]

## _TODO_ diameter or radii for unfolding (see Wicksel)?
ret <- unfold(S$sp,nclass=25)
 
## Point process intensity
cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")
 
## original diameters
r3d <- unlist(lapply(S,function(x) 2.0*x$r))
rest3d <- unlist(lapply(2:(length(ret$breaks)),
            function(i) rep(ret$breaks[i],sum(ret$N_V[i-1]))))
 
op <- par(mfrow = c(1, 2))
hist(r3d[r3d<=max(ret$breaks)], breaks=ret$breaks, main="Radius 3d",
     freq=FALSE, col="gray",xlab="r")
hist(rest3d, breaks=ret$breaks,main="Radius estimated",
     freq=FALSE, col="gray", xlab="r")
par(op)



