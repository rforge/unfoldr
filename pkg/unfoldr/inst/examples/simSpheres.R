# Comment: Simulate Poisson sphere system and intersect 
# 
# Author: bahama
###############################################################################

library(rgl)
library(plotrix)
library(unfoldr)

spheres <- function(spheres, box=NULL, draw.box=FALSE, draw.axis=FALSE, ...) {
	xyz <- apply(sapply(spheres, "[[", "center"),1,function(x) x)
	sizes <- unlist(lapply(spheres,function(x) x$r))
	spheres3d(xyz,radius=sizes,...)
	
	if(draw.axis) {
		axes3d(c('x','y','z'), pos=c(1,1,1), tick=FALSE)
		#title3d('','','x','y','z')		
	}
	axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
			expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
	
	if(draw.box) {
		x <- box$xrange[2]
		y <- box$yrange[2]
		z <- box$zrange[2]	
		c3d.origin <- rgl::translate3d(rgl::scale3d(rgl::cube3d(col="gray", alpha=0.1),x/2,y/2,z/2),x/2,y/2,z/2)
		rgl::shade3d(c3d.origin)
	}
	
}

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

#################################################
## beta distribution for radii
#################################################

lam <- 50
## parameter beta distribution (radii)
theta <- list("size"=list("shape1"=1,"shape2"=10))

# simulation bounding box
box <- list("xrange"=c(-1,4),"yrange"=c(-1.5,3.5),"zrange"=c(0,2))

## simulate and return full spheres system
## intersect with XZ plane and return full list of intersection profiles

## section plane xy
S <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",n=c(0,0,1),dz=0,pl=101)

# check resulting distribution
length(S$S)
summary(sapply(S$S,"[[","r"))
theta$size[[1]]/(theta$size[[1]]+theta$size[[2]])			# mean

## interior spheres:
## the ones which intersect one of the lateral planes (without top/bottom planes)
## showing spheres with color intersect 
notIn <- sapply(S$S,function(x) !attr(x,"interior"))
spheres(S$S[notIn],box,TRUE,color=col)

# not intersecting
In <- sapply(S$S,function(x) attr(x,"interior"))
spheres(S$S[In],box,color="gray")

## ful sphere system
#open3d()
#spheres(S$S,box,TRUE,color=col)
#planes3d(0,0,1,0,col="darkgray",alpha=1)

## draw intersections
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
spheres(S$S[id],box,TRUE,color=col)
planes3d(0,0,1,0,col="darkgray",alpha=1)

XYr <- t(sapply(sp,function(s) cbind(s$center[1],s$center[2],s$r)))
# centers
x <- XYr[,1]
y <- XYr[,2]
r <- XYr[,3]
xlim <- c(-1,4)
ylim <- c(-1.5,3.5) 

dev.new()
plot(x,y,type="n",xaxs="i", yaxs="i", xlab="x",ylab="y",xlim=xlim,ylim=ylim)
for(i in 1:nrow(XYr))
 draw.circle(x[i],y[i],r[i],nv=100,border=NULL,col="black")

## digitize inersections
sp <- S$sp
## `win` can also be omitted if not different
## from original itersection window (derived from the box)
win <- attr(sp,"win")
W <- digitizeProfiles(sp, delta=0.01, win=win)
dim(W)
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

######################################################
## exact simulation of spheres with log normal radii
######################################################

lam <- 100
## parameter rlnorm distribution (radii)
theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5))
# simulation bounding box
box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
## simulate only 3D system
S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="sphere",
		intersect="original", perfect=TRUE,pl=101)
## show
spheres(S[1:2000],box,TRUE,TRUE,color=col)

## check!
mean(log(sapply(S,"[[","r")))
sd(log(sapply(S,"[[","r")))


#######################################################
## planar section
#######################################################

# planar section of exact simulated `rlnorm` sphere system:
# returns diameters for those section profiles
# which have their centers inside the intersection window
sp <- planarSection(S,d=2.5,intern=TRUE,pl=1)
hist(sp)
summary(sp)
mean(log(sp/2))
sd(log(sp/2))

#######################################################
## unfolding
#######################################################
ret <- unfold(sp,nclass=25)

## Point process intensity
cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")

## original diameters
r3d <- unlist(lapply(S,function(x) 2*x$r))
rest3d <- unlist(lapply(2:(length(ret$breaks)),
				function(i) rep(ret$breaks[i],sum(ret$N_V[i-1]))))

op <- par(mfrow = c(1, 2))
hist(r3d[r3d<=max(ret$breaks)], breaks=ret$breaks, main="Radius 3d",
		freq=FALSE, col="gray",xlab="r")
hist(rest3d, breaks=ret$breaks,main="Radius estimated",
		freq=FALSE, col="gray", xlab="r")
par(op)


