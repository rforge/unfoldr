# Comment: Simulate Poisson sphere system and intersect 
# 
# Author: bahama
###############################################################################

library(rgl)
library(plotrix)
library(unfoldr)

# drawing function

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

## intensity
lam <- 500

# simulation bounding box
box <- list("xrange"=c(0,2),"yrange"=c(0,2),"zrange"=c(0,5))

## show how to call
#head(simPoissonSystem)

.tA <- function(p,ctr,a,b,phi) {
	(p[1]-ctr[1] + sin(phi)/cos(phi)*(p[2]-ctr[2]))/(a*cos(phi)+a*((sin(phi))^2)/cos(phi))
}

.PointOnEllipse <- function(ctr,a,b,phi,tt) {
	pxy <- numeric(2)
	pxy[1] = ctr[1] + a*cos(tt)*cos(phi)-b*sin(tt)*sin(phi)
	pxy[2] = ctr[2] + a*cos(tt)*sin(phi)+b*sin(tt)*cos(phi)
	return (pxy)
}

draw.EllipseSector <- function(x, ctr , nv=100,...) {
	
	psi0_1 <- ifelse(x$pS<0,-x$psi[1],x$psi[1])
	psi1_1 <- ifelse(x$pS<0, x$psi[1],2*pi-x$psi[1])
	
	psi0_2 <- ifelse(x$pS>0,-x$psi[2],x$psi[2])
	psi1_2 <- ifelse(x$pS>0, x$psi[2],2*pi-x$psi[2])
	
	tt <- sort(c(psi0_1,psi1_1,psi0_2,psi1_2))
	
	pp1 <- rbind(t(sapply(seq(tt[2],tt[3],length=nv), function(xi) .PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,xi))))
	pp2 <- rbind(t(sapply(seq(tt[4],2*pi-abs(tt[1]),length=nv), function(xi) .PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,xi))))
	
	pp <- rbind(.PointOnEllipse(ctr,x$ab[1],x$ab[2],x$phi,tt[1]), pp1 ,pp2)
	polygon(pp[,1], pp[,2],...)
	
}

draw.segments <- function(E, cyltype=0, x=c(0,1), y=x, angle=0, normal=c(0,0,1),axes=TRUE, ...)
{
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
					draw.EllipseSector(x,x$center[c(i,j)], ...)
					
					## fist end cap
					t <- x$psi[2]
					## get angle on circle - not ellipse!
					p <- .PointOnEllipse(x$center[c(i,j)],x$ab[1],x$ab[2],x$phi,t)
					t0 <- acos(.tA(p,x$mPoint1[c(i,j)],x$rcaps[2],x$rcaps[2],x$phi))
					psi0 <- ifelse(x$pS>0,-t0,t0)
					psi1 <- ifelse(x$pS>0, t0,2*pi-t0)
					draw.ellipse(x$mPoint1[i],x$mPoint1[j], x$rcaps[2], x$rcaps[2],angle=x$phi,	segment=c(psi0,psi1), deg=FALSE,...)
					
					## second end cap
					t <- x$psi[1]
					## get angle on circle - not ellipse!
					p <- .PointOnEllipse(x$center[c(i,j)],x$ab[1],x$ab[2],x$phi,t)
					t0 <- acos(.tA(p,x$mPoint0[c(i,j)],x$rcaps[1],x$rcaps[1],x$phi))
					psi0 <- ifelse(x$pS<0,-t0,t0)
					psi1 <- ifelse(x$pS<0, t0,2*pi-t0)
					draw.ellipse(x$mPoint0[i],x$mPoint0[j], x$rcaps[1], x$rcaps[1],angle=x$phi,segment=c(psi0,psi1), deg=FALSE,...)
						
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
					#phi <- cartesian2sphere(xx$major)
					
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
						
						#pxy <- .PointOnCircle(x$mPoint0,x$rcaps[1],x$phi,ifelse(x$pS<0,0,pi))
						#points(pxy[1],pxy[2],col="red")
						#segments(x$ipt0[1],x$ipt0[2], pxy[1],pxy[2],col="red")
						
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
## simulate bivariate
#################################################################

## no `shape` required
theta <- list("size"=list("mx"=-1.45,"my"=-2.0, "sdx"=0.15,"sdy"=0.25,"rho"=0.0),
			  "orientation"=list("kappa"=5))

S <- simPoissonSystem(theta,lam,size="rbinorm",box=box,type="cylinders",
		intersect="full", n=c(0,0,1), "orientation"="rbetaiso", dz=2.5,
		perfect=TRUE,intern=FALSE,delta=0.005, pl=101)

#open3d()
#cylinders3d(S$S[1:1000], draw.box=TRUE, box=box, col=col)
#planes3d(0,0,-1,2.5,col="black",alpha=1)

## 3D intersected objects
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
cylinders3d(S$S[id], draw.box=TRUE, box=box, col=col)
planes3d(0,0,-1,2.5,col="black",alpha=1)

draw.segments(sp,x=c(0,2),bg="gray",normal=c(0,0,1),xlab="[mm]",xaxs="i",yaxs="i",
		ylab="[mm]",col="black",axes=TRUE,cex.lab=1.8,cex=1.8,cex.axis=1.8)

## digitized
W <- S$W
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

##################################################################################
# general intersections (should be same as above)
##################################################################################

S <- S$S
SP <- intersectSystem(S, 2.5, n=c(0,0,1), intern=FALSE, pl=1)
id <- sapply(SP,"[[","id") 
length(id)

open3d()
cylinders3d(S[id], draw.box=TRUE, box=box, col=col)
planes3d(0,0,-1,2.5,col="black",alpha=1)

draw.segments(SP,x=c(0,2),bg="gray",normal=c(0,0,1),xlab="[mm]",xaxs="i",yaxs="i",
		ylab="[mm]",col="black",axes=TRUE,cex.lab=1.8,cex=1.8,cex.axis=1.8)


## digitize
W2 <- digitizeProfiles(SP, delta=0.005)
dev.new()
image(1:nrow(W2),1:ncol(W2),W2,col=gray(1:0))


#################################################################
## Update intersection: find objects which intersect bounding box
#################################################################

idx <- updateIntersections(S)
sum(!idx)							# objects intersecting
id <- which( idx != 1)	

# show in 3D
open3d()
spheroids3d(S[id], FALSE, TRUE, box=box, col=col)

#################################################################
## user-defined simulation function
#################################################################


# no perfect simualtion here for 'rmulti'
# multivariate size distribution,
# independent orientation distribution 
rmulti <- function(m,s,kappa) {	
	# directional distribution
	# (implemented `rbetaiso` distribution)
	rbetaiso <- function(kappa) {
		phi <- runif(1,0,1)*2*pi
		q <- runif(1,0,1)
		theta=acos((1-2*q)/sqrt(kappa*kappa-(1-2*q)*(1-2*q)*(kappa*kappa-1)))
		list("u"=c(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)),
				"theta"=theta,"phi"=phi)					
	}
	
	dir <- rbetaiso(kappa)
	# log normal semi-major/semi-minor lengths
	M <- chol(s, pivot = TRUE)
	M <- M[,order(attr(M, "pivot"))]
	x <- exp(matrix(m,nrow=1) + matrix(rnorm(ncol(s)), nrow = 1, byrow = TRUE) %*% M)
	a <- min(x)
	b <- max(x)
	# the following elements are obligatory as a
	# return value for user-defined spheroid simulations
	list("a"=a,"b"=b,"c"=a,"u"=dir$u,"shape"=a/b,"theta"=dir$theta, "phi"=dir$phi)	
}


sigma <- matrix(c(0.1,0.1,0.1,0.25), ncol=2)
theta <- list("m"=c(-3.0,-2.0),"s"=sigma,"kappa"=0.5)

S <- simPoissonSystem(theta,lam,rjoint=rmulti,box=box,type="prolate",
		intersect="full",n=c(0,0,1), mu=c(0,0,1), dz=2.5, pl=101)

# in 3D
sp <- S$sp # sections
id <- sapply(sp,"[[","id") 
open3d()
spheroids3d(S$S[id], FALSE, TRUE, box=box, col=col)
planes3d(0,0,-1,2.5,col="black",alpha=1)

# in 2D
dev.new()
Es <- drawEllipses(sp, x=box$xrange, border="black",xlab="[mm]", ylab="[mm]",
		rot=0, bg="gray",col=col,	cex.lab=1.8,cex=1.8,cex.axis=1.8,nv=1000)

# digitized
W <- S$W
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))
