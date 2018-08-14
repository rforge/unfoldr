\dontrun{
## Comment: Simulate a Poisson spheroid system,
## 			intersect, discretize and display results
	
library(rgl)
library(plotrix)

drawEllipses <- function(E, x=c(0,1), y=x, xlab="x",ylab="y", bg="gray", angle=0, ...) {
	Es <- sapply(E,
			function(x) {
				c("id"=x$id,
				  "x"=x$center[1],
				  "y"=x$center[2],
				  "phi"=x$phi,
				  "a"=x$ab[1],
				  "b"=x$ab[2])
			})
	
	plot(x,y, type="n",xlab=xlab,ylab=ylab)#,xaxs="i", yaxs="i")
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bg)
	
	if(angle>0) {
		M <- matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2 )
		XY <- matrix( c(Es["x",],Es["y",]), nc=2)  %*% M
		draw.ellipse(-XY[,1],XY[,2], Es["a",], Es["b",],
			angle=Es["phi",]+angle, deg=FALSE,...)
	} else  {
		draw.ellipse(Es["x",], Es["y",], Es["a",], Es["b",],
			angle=Es["phi",], deg=FALSE,...)
	}
	
	ret <- apply(Es,2,function(x) as.list(x))
	return ( ret )
}

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

## intensity
lam <- 100

# simulation bounding box
box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))

# show how to call
head(simPoissonSystem)

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
spheroids3d(S[1:1000], FALSE, TRUE, box=box, col=col)

###################################################################
## simulate bivariate size-shape distribution for prolate spheroids
###################################################################

## no `shape` required, isotropic orientation
theta <- list("size"=list("mx"=-2.5,"my"=0.5, "sdx"=0.35,"sdy"=0.25,"rho"=0.15),
		"orientation"=list("kappa"=1))

S <- simPoissonSystem(theta,lam,size="rbinorm",box=box,type="prolate",
		intersect="full", ,n=c(0,0,1), mu=c(0,0,1),
		"orientation"="rbetaiso", dz=2.5,perfect=TRUE,pl=101)

## 3D intersected objects
sp <- S$sp
id <- sapply(sp,"[[","id") 
open3d()
spheroids3d(S$S[id], FALSE, TRUE, box=box, col=col)
planes3d(0,0,-1,2.5,col="black",alpha=1)

## angle always w.r.t to x axis after simulation
phi <- sapply(sp,"[[","phi")
summary(phi)

## check rotation matrix ellipses
#E <- sp[[10]]
#V <- eigen(E$A)
#1/sqrt(V$values)
#
#A <- diag(c(1/E$ab[1]^2,1/E$ab[2]^2))
#B <- cbind(E$major,E$minor)
#t(B) %*% A %*% B
#E$A

dev.new()
Es <- drawEllipses(sp, x=box$xrange, y=box$yrange, border="black",xlab="[mm]", ylab="[mm]",
		bg="gray",col=col,	cex.lab=1.8,cex=1.8,cex.axis=1.8,nv=1000)

## digitized image
W <- S$W
dev.new()
image(1:nrow(W),1:ncol(W),W,col=gray(1:0))

##################################################################################
# general intersections (should be same as above)
##################################################################################

S <- S$S
SP <- intersectSystem(S, 2.5, n=c(0,0,1), intern=FALSE, pl=1)

## check compare 
str(sp[[100]])
str(SP[[100]])

#################################################################
## Vertical section
#################################################################

# vertical intersection
dz <- 2.5		# distance to origin of box [0,5]^3
n <- c(0,1,0)	# normal in y direction (xz plane)

theta$orientation$kappa <- 10	# random planar in xy plane

S <- simPoissonSystem(theta,lam,size="rbinorm",box=box,
		type="prolate", intersect="full", ,n=n, mu=c(0,0,1),
		"orientation"="rbetaiso", dz=dz, perfect=TRUE, intern=TRUE, pl=101)

sp <- S$sp # sections
id <- sapply(sp,"[[","id") 
open3d()
spheroids3d(S$S[id], FALSE, TRUE, box=box, col=col)
planes3d(0,-1,0,2.5,col="black",alpha=1)

# check angle
phi <- sapply(sp,"[[","phi")
summary(phi)

dev.new()
Es <- drawEllipses(sp, x=box$xrange, y=box$yrange, border="black",xlab="[mm]", ylab="[mm]",
		 bg="gray",col=col,	cex.lab=1.8,cex=1.8,cex.axis=1.8,nv=1000)
 
# intersect 3D system
Sp <- S$S   # spheroids
spv <- verticalSection(Sp,d=dz,n=n,intern=TRUE)

# angle in the intersecting plane always
# w.r.t. to mu=(0,0,1) used for unfolding
summary(spv$alpha) 

#################################################################
## Update intersection: find objects which intersect bounding box
#################################################################

idx <- updateIntersections(Sp)
sum(!idx)							# objects intersecting
id <- which( idx != 1)	

# show in 3D
open3d()
spheroids3d(Sp[id], FALSE, TRUE, box=box, col=col)

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
}
