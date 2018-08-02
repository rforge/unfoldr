# Comment: Simulate Poisson sphere system and intersect 
# 
# Author: bahama
###############################################################################

library(rgl)
library(plotrix)
library(unfoldr)

drawEllipses <- function(E, x=c(0,1), y=x, xlab="x",ylab="y", bg="gray", angle=0, rot=0, ...) {
	Es <- sapply(E,
			function(x) {
				c("id"=x$id,
				  "x"=x$center[1],
				  "y"=x$center[2],
				  "phi"=x$phi + rot,		# add pi/2 to conform with `u` direction
				  "a"=x$ab[1],
				  "b"=x$ab[2])
			})
	
	plot(x,y, type="n",xlab=xlab,ylab=ylab)#,xaxs="i", yaxs="i")
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bg)
	
	if(angle>0) {
		M <- matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2 )
		XY <- matrix( c(Es["x",],Es["y",]), nc=2)  %*% M
		draw.ellipse(-XY[,1],XY[,2], Es["a",], Es["b",], angle=Es["phi",]+angle, deg=FALSE,...)
	} else  {
		draw.ellipse(Es["x",], Es["y",], Es["a",], Es["b",], angle=Es["phi",], deg=FALSE,...)
	}
	
	ret <- apply(Es,2,function(x) as.list(x))
	return ( ret )
}

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

## intensity
lam <- 100

# simulation bounding box
box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))


#################################################################
## simulate bivariate
#################################################################

## no `shape` required
theta <- list("size"=list("mx"=-2.5,"my"=0.5, "sdx"=0.35,"sdy"=0.25,"rho"=0.15),
		"orientation"=list("kappa"=10))

head(simPoissonSystem)
S <- simPoissonSystem(theta,lam,size="rbinorm",box=box,type="prolate",
		intersect="full", ,n=c(0,0,1), mu=c(0,0,1),
		"orientation"="rbetaiso", dz=2.5,perfect=TRUE,pl=101)

#theta <- sapply(S$S,function(x) x$angles[1])

## check!
#ACB <- t(sapply(S$S,"[[","acb"))
#mean(log(ACB[,3])) # mx
#sd(log(ACB[,3]))   # sdx

## show
#open3d()
#spheroids3d(S$S[1:2000], FALSE, TRUE, box=box, col=col)

## 3D intersected objects
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
spheroids3d(S$S[id], FALSE, TRUE, box=box, col=col)
planes3d(0,-1,0,2.5,col="black",alpha=1)
#planes3d(-1,0,0,2.5,col="black",alpha=1)
planes3d(0,0,-1,2.5,col="black",alpha=1)

(phi <- sapply(sp,"[[","phi"))
summary(phi)

## as 2D intersecions
#sp2 <- sp
#for(i in 1:length(sp))
#	sp2[[i]]$phi <- .getAngle(sp2[[i]]$phi)

#E <- sp[[10]]
#V <- eigen(E$A)
#1/sqrt(V$values)
#
#A <- diag(c(1/E$ab[1]^2,1/E$ab[2]^2))
#B <- cbind(E$major,E$minor)
#t(B) %*% A %*% B

dev.new()
Es <- drawEllipses(sp, x=box$xrange, border="black",xlab="[mm]", ylab="[mm]",
		rot=0, bg="gray",col=col,	cex.lab=1.8,cex=1.8,cex.axis=1.8,nv=1000)


## digitized
W <- S$W
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

####################################################################
## `rlnorm` distribution 
####################################################################

# log normal size, constant shape, isotropic orientation (rbetaiso) 
theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5),
			  "shape"=list("s"=0.5),
			  "orientation"=list("kappa"=1))

## simulate and return full spheres system
## intersect with XZ plane and return full list of intersection profiles
S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="oblate",
		intersect="original",mu=c(0,1,0),n=c(0,0,1),dz=0,perfect=TRUE,pl=101)

## show some objects to get an impression
open3d()
spheroids3d(S[1:2000], FALSE, TRUE, box=box, col=col)


#################################################################


## general intersections (same as above)
## but only return interior section profiles
S <- S$S
SP <- intersectSystem(S, 2.5, n=c(0,1,0), intern=FALSE, pl=1)

## check compare 
str(sp[[100]])
str(SP[[100]])

## show in 3D (section profiles at box front)
id <- sapply(SP,"[[","id") 
open3d()
spheroids3d(S[id], FALSE, TRUE, box=box, col=col)
planes3d(0,-1,0,2.5,col="darkgray",alpha=1)
#drawSpheroidIntersection(SP,n=c(0,1,0),np=20)

#phi2 <- sapply(sp,"[[","phi")
#phi <- sapply(SP,"[[","phi")

dev.new()
EsIntern <- drawEllipses(SP, x=box$xrange, border="black",xlab="[mm]", ylab="[mm]",
				rot=0,bg="gray",col=col,cex.lab=1.8,cex=1.8,cex.axis=1.8,nv=1000)


## digitize 2D sections
W <- digitizeProfiles(SP, delta=0.01)
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

dev.new()
image(t(W))

## TODO:

#################################################################
## Vertical section
#################################################################


#################################################################
## Unfolding
#################################################################



#################################################################
## Update intersection: find objects which intersect bounding box
#################################################################



#################################################################
## user-defined simulation function
#################################################################