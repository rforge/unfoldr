\dontrun{
## Comment: Simulate a Poisson spherocylinder system,
## 			intersect, discretize and display results

library(rgl)
library(plotrix)

## intensity
lam <- 500

## simulation bounding box
box <- list("xrange"=c(0,2),"yrange"=c(0,2),"zrange"=c(0,5))
col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")

## 2D intersection of spherocylinders
draw.segments <- function(E, cyltype=0, x=c(0,1), y=x, angle=0, normal=c(0,0,1),axes=TRUE, ...)
{
	.tA <- function(p,ctr,a,b,phi) {
		(p[1]-ctr[1] + sin(phi)/cos(phi)*(p[2]-ctr[2]))/(a*cos(phi)+a*((sin(phi))^2)/cos(phi))
	}
	
	.PointOnEllipse <- function(ctr,a,b,phi,tt) {
		pxy <- numeric(2)
		pxy[1] = ctr[1] + a*cos(tt)*cos(phi)-b*sin(tt)*sin(phi)
		pxy[2] = ctr[2] + a*cos(tt)*sin(phi)+b*sin(tt)*cos(phi)
		return (pxy)
	}
	
	.EllipseSector <- function(x, ctr , nv=100,...)
	{
		psi0_1 <- ifelse(x$pS<0,-x$psi[1],x$psi[1])
		psi1_1 <- ifelse(x$pS<0, x$psi[1],2*pi-x$psi[1])
		
		psi0_2 <- ifelse(x$pS>0,-x$psi[2],x$psi[2])
		psi1_2 <- ifelse(x$pS>0, x$psi[2],2*pi-x$psi[2])
		
		tt <- sort(c(psi0_1,psi1_1,psi0_2,psi1_2))
		
		pp1 <- rbind(t(sapply(seq(tt[2],tt[3],length=nv),
				function(xi) .PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,xi))))
		pp2 <- rbind(t(sapply(seq(tt[4],2*pi-abs(tt[1]),length=nv),
				function(xi) .PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,xi))))
		
		pp <- rbind(.PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,tt[1]), pp1 ,pp2)
		polygon(pp[,1], pp[,2],...)	
	}
	
	
	argg <- c(as.list(environment()), list(...))
	mycol <- "black"
	if(!is.null(argg$col) && length(argg$col) > 0)
		mycol<-argg$col
	
	plot(x,y,type="n",axes=axes,...)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray")
	
	## plane indices
	i <- j <- 0
	switch(which(normal==1, arr.ind=TRUE),{i<-2;j<-3}, {i<-1;j<-3},{i<-1;j<-2})
	
	idType <- unlist(lapply(E, function(x) if( x$type==9) TRUE else FALSE) )
	if(length(idType)>0 && any(idType)) {
		## draw caps
		X <- E[idType]
		lapply(1:length(X), function(k) {
					x <- X[[k]]
					.EllipseSector(x,x$center[c(i,j)], ...)
					
					## fist end cap
					t <- x$psi[2]
					## get angle on circle - not ellipse!
					p <- .PointOnEllipse(x$center[c(i,j)],x$ab[1],x$ab[2],x$phi,t)
					t0 <- acos(.tA(p,x$mPoint1[c(i,j)],x$rcaps[2],x$rcaps[2],x$phi))
					psi0 <- ifelse(x$pS>0,-t0,t0)
					psi1 <- ifelse(x$pS>0, t0,2*pi-t0)
					draw.ellipse(x$mPoint1[i],x$mPoint1[j], x$rcaps[2], x$rcaps[2],
						angle=x$phi,segment=c(psi0,psi1), deg=FALSE,...)
					
					## second end cap
					t <- x$psi[1]
					## get angle on circle - not ellipse!
					p <- .PointOnEllipse(x$center[c(i,j)],x$ab[1],x$ab[2],x$phi,t)
					t0 <- acos(.tA(p,x$mPoint0[c(i,j)],x$rcaps[1],x$rcaps[1],x$phi))
					psi0 <- ifelse(x$pS<0,-t0,t0)
					psi1 <- ifelse(x$pS<0, t0,2*pi-t0)
					
					draw.ellipse(x$mPoint0[i],x$mPoint0[j], x$rcaps[1], x$rcaps[1],
						angle=x$phi,segment=c(psi0,psi1), deg=FALSE,...)						
				})
	}
	
	## Iterate over all but Ellipse segments
	idType <- unlist(lapply(E, function(x) if( x$type!=9) TRUE else FALSE) )
	if(length(idType)>0 && any(idType))
	{
		Es <- sapply(E[idType], function(x) {
					psi0 <- 0
					psi1 <- 2*pi
					
					if(x$type==5){
						## Circle as intersection
						ret <- c("type"=x$type,"id"=x$id, "x"=x$mPoint0[i],"y"=x$mPoint0[j],"phi"=0,
								"a"=x$radius, "b"=x$radius,"psi0"=psi0, "psi1"=psi1)
					} else	if(x$type==6) {
						## Cap
						ret <- c("type"=x$type,"id"=x$id, "x"=x$mPoint0[i],"y"=x$mPoint0[j],"phi"=0,
								"a"=x$rcaps[1], "b"=x$rcaps[1],"psi0"=psi0, "psi1"=psi1)
					} else if(x$type==7) {
						## Ellipse
						ret <- c("type"=x$type,"id"=x$id, "x"=x$center[i],"y"=x$center[j],"phi"=x$phi,
								"a"=x$ab[1], "b"=x$ab[2],"psi0"=psi0, "psi1"=psi1)
					} else if(x$type==8) {
						## Ellipse arcs
						psi0 <- ifelse(x$pS>0,-x$psi[1],x$psi[1])
						psi1 <- ifelse(x$pS>0,x$psi[1],2*pi-x$psi[1])
						ret <- c("type"=x$type,"id"=x$id, "x"=x$center[i],"y"=x$center[j],"phi"=x$phi,
								"a"=x$ab[1], "b"=x$ab[2],"psi0"=psi0, "psi1"=psi1)
					}
					ret
				}
		)
		
		draw.ellipse(Es["x",], Es["y",], Es["a",], Es["b",], angle=Es["phi",],
				segment=cbind(Es["psi0",],Es["psi1",]), deg=FALSE,...)
		
	}
	
	## Ellipse arcs
	idType <- unlist(lapply(E, function(x) if( x$type==8) TRUE else FALSE) )
	if(length(idType)>0 && any(idType)) {
		X <- E[idType]
		lapply(1:length(X), function(k) {
					x <- X[[k]]
					if(x$rcaps[1]>0) {
						# plot circle
						# get point on ellipse where to cut off and
						# get the cut off angle for the circle
						t <- x$psi[1]
						p <- .PointOnEllipse(x$center[c(i,j)],x$ab[1],x$ab[2],x$phi,t)
						t0 <- acos(.tA(p,x$mPoint0[c(i,j)],x$rcaps[1],x$rcaps[1],x$phi))
						psi0 <- ifelse(x$pS<0,-t0,t0)
						psi1 <- ifelse(x$pS<0, t0,2*pi-t0)
						
						draw.ellipse(x$mPoint0[i],x$mPoint0[j],x$rcaps[1],x$rcaps[1],
								angle=x$phi,segment=c(psi0,psi1), deg=FALSE,...)
		
					}
				})
	}	
}

####################################################################
## `rlnorm` distribution 
####################################################################

# log normal size, constant shape, isotropic orientation (rbetaiso) 
theta <- list("size"=list("meanlog"=-1.45,"sdlog"=0.15),
		"shape"=list("s"=0.25), "orientation"=list("kappa"=1))

## simulate and return full spheres system
## intersect with XZ plane and return full list of intersection profiles
S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="cylinders",
		intersect="original",n=c(0,0,1),dz=2.5,pl=101)

## show some objects to get an impression
open3d()
cylinders3d(S[1:1000], draw.box=TRUE, box=box, col=col)


#################################################################
## simulate bivariate size/shape distribution
#################################################################

## no extra `shape` factor required, see documentation
theta <- list("size"=list("mx"=-1.45,"my"=-2.0, "sdx"=0.15,"sdy"=0.25,"rho"=0.0),
			  "orientation"=list("kappa"=5))

S <- simPoissonSystem(theta,lam,size="rbinorm",box=box,type="cylinders",
		intersect="full", n=c(0,1,0), "orientation"="rbetaiso", dz=1.5,
		perfect=TRUE, intern=TRUE, delta=0.005, pl=101)

#open3d()
#cylinders3d(S$S[1:1000], draw.box=TRUE, box=box, col=col)
#planes3d(0,0,-1,2.5,col="black",alpha=1)

## 3D intersected objects
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
cylinders3d(S$S[id], draw.box=TRUE, box=box, col=col)
planes3d(0,-1,0,1.5,col="black",alpha=1)

n <- attr(sp,"plane")	# normal vector of intersecting plane
win <- attr(sp,"win")	# intersection window

dev.new()
draw.segments(sp,x=win[[1]],y=win[[2]],bg="gray",normal=n,xlab="[mm]",xaxs="i",yaxs="i",
		ylab="[mm]",col="black",axes=TRUE,cex.lab=1.8,cex=1.8,cex.axis=1.8)

## digitized
W <- S$W
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

##################################################################################
# general intersections (should be same as above)
##################################################################################

Sp <- S$S
SP <- intersectSystem(Sp, 1.5, n=c(0,1,0), intern=TRUE, pl=1)
id <- sapply(SP,"[[","id") 

open3d()
cylinders3d(Sp[id], draw.box=TRUE, box=box, col=col)
planes3d(0,-1,0,1.5,col="black",alpha=1)

n <- attr(SP,"plane")	# normal vector of intersecting plane
win <- attr(SP,"win")	# intersection window

dev.new()
draw.segments(SP,x=win[[1]],y=win[[2]],bg="gray",normal=n,xlab="[mm]",xaxs="i",yaxs="i",
		ylab="[mm]",col="black",axes=TRUE,cex.lab=1.8,cex=1.8,cex.axis=1.8)


## digitize
W2 <- digitizeProfiles(SP, delta=0.005)
dev.new()
image(1:nrow(W2),1:ncol(W2),W2,col=gray(1:0))

#################################################################
## Update intersection: find objects which intersect bounding box
#################################################################

idx <- updateIntersections(Sp)
sum(!idx)							# number of objects intersecting bounding box
id <- which( idx != 1)	

## show in 3D
# open3d()
# cylinders3d(Sp[id], draw.box=TRUE, box=box, col=col)


#################################################################
## example of a user-defined distributino of size and shape
#################################################################

# no perfect simualtion here:
# log normal length of cylinder, beta distributed shape factor,
# independent orientation distribution 
rmulti <- function(m,s) {	
	# directional distribution	
	phi <- runif(1)*2*pi
	theta <- acos(2.0*runif(1)-1.0)
	dir <- list("u"=c(cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)),
			    "theta"=theta,"phi"=phi)
	len <- rlnorm(1,m,s)
	shape <- rbeta(1,1,10)
	r <- len*shape
	
	# the following elements are obligatory as return values within a list
	list("h"=len-2*r,"r"=r,"u"=dir$u,"shape"=shape,"theta"=dir$theta, "phi"=dir$phi)	
}

# parameter for size distribution
theta <- list("m"=-2.25,"s"=0.35)
# simulate system
S <- simPoissonSystem(theta,lam,rjoint=rmulti,box=box,type="cylinders",
		intersect="full",n=c(0,0,1), mu=c(0,0,1), dz=2.5, pl=101)

# in 3D
sp <- S$sp # sections
id <- sapply(sp,"[[","id") # indices of intersected objects 
# show
open3d()
cylinders3d(S$S[id], draw.box=TRUE, box=box, col=col)
planes3d(0,0,-1,2.5,col="black",alpha=1)

# 2D intersection profiles
n <- attr(sp,"plane")	# normal vector of intersecting plane
win <- attr(sp,"win")	# intersection window

dev.new()
draw.segments(sp,x=win[[1]],y=win[[2]],bg="gray",normal=n,xlab="[mm]",
	xaxs="i",yaxs="i",ylab="[mm]",col="black",axes=TRUE,cex.lab=1.8,cex=1.8,cex.axis=1.8)

## digitize objects
W <- digitizeProfiles(sp, delta=0.001)
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))
}
