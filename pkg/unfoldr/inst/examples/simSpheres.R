# Comment: Simulate Poisson sphere system and intersect 
# 
# Author: bahama
###############################################################################

library(rgl)
library(plotrix)
library(unfoldr)

spheres <- function(spheres, box=NULL, draw.axis=FALSE, ...) {
	xyz <- apply(sapply(spheres, "[[", "center"),1,function(x) x)
	sizes <- unlist(lapply(spheres,function(x) x$r))
	spheres3d(xyz,radius=sizes,...)
	
	if(draw.axis) {
		axes3d(c('x','y','z'), pos=c(1,1,1), tick=FALSE)
		#title3d('','','x','y','z')		
	}
	axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
			expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
	
}

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

## beta distribution for radii
lam <- 50
## parameter beta distribution (radii)
theta <- list("size"=list("shape1"=1,"shape2"=10))

# simulation bounding box
box <- list("xrange"=c(-1,4),"yrange"=c(-1.5,3.5),"zrange"=c(0,2))

## simulate and return full spheres system
## intersect with XZ plane and return full list of intersection profiles
S <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",n=c(0,1,0),dz=-1.5,pl=101)


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
open3d()
spheres(S$S,box,TRUE,color=col)
planes3d(0,-1,0,-1.5,col="darkgray",alpha=1)

## draw intersections
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
spheres(S$S[id],box,TRUE,color=col)
planes3d(0,-1,0,-1.5,col="darkgray",alpha=1)

XYr <- t(sapply(sp,function(s) cbind(s$center[1],s$center[3],s$r)))
# centers
x <- XYr[,1]
y <- XYr[,2]
r <- XYr[,3]
xlim <- c(-1,4)
ylim <- c(0,2)


dev.new()
plot(x,y,type="n",xaxs="i", yaxs="i", xlab="x",ylab="y",xlim=xlim,ylim=ylim)
for(i in 1:nrow(XYr))
 draw.circle(x[i],y[i],r[i],nv=100,border=NULL,col="black")

## digitize inersections
sp <- S$sp
win <- attr(sp,"win")

W <- digitizeProfiles(sp, delta=0.01, win=win)
dim(W)

dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))
