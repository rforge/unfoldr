###############################################################################
# Author:  M. Baaske
# Date:	   2018/06/15	
# File:    spheroid.R: 
# 
# Comment: simulation, intersections and visualization of spheroid systems
# 
###############################################################################

## angle in the section plane relative to z axis in 3D
.getAngle <- function(phi) {
	## this is slower!
	# return (abs(asin(sin(s))))
	
	if(phi<=pi/2) { phi }
	else {
		if( phi <= pi) pi-phi
		else if(phi<1.5*pi) phi %% pi
		else 2*pi-phi
	}
}


#' Check intersections with simultion box
#'
#' Check itersections of either spheres, spheroids or cylinders with the bounding (simulation) box
#'
#' For a given list of spheres, spheroids or cylinders the function tests whether an object intersects
#' the bounding simulation box. The vector returned is of length equal to \code{S} and has entries either \code{1}
#' for an object which is intersected by one of the bounding lateral planes of the simulation box or otherwise \code{0}.
#'
#' @param S 	list of spheres, spheroids or cylinders, see \code{\link{simPoissonSystem}}
#'  
#' @return 		binary integer vector of length equal to the length of \code{S}.
#' 			    THe coding is \code{0} for not interior and \code{1} for interior 
#'  
#' @author M. Baaske
#' 
#' @rdname updateIntersections
#' @export
updateIntersections <- function(S) {
	.Call(C_UpdateIntersections,as.character(substitute(S)), .GlobalEnv)
}

#' Construct section profiles of spheroids
#'
#' Get a storing structure for section profiles of spheroids
#'
#' The function aggregates the necessary data for the trivariate unfolding procedure either for \code{type}
#' "\code{prolate}" or "\code{oblate}" spheroids. Argument \code{size} is a numeric matrix which contains the
#' axes lengths (first column corresponds to major semi-axis, second one to minor semi-axis). The \code{angle}
#' is the orientation angle between the major axis and the vertical axis direction in the vertical intersection plane.
#' If the angles range within \eqn{[0,2\pi]} these are transformed to \eqn{[0,\pi/2]}. The function returns a
#' list which consists of either longer or shorter axis \code{A} of section profiles corresponding to the type
#' of spheroids which are intended to be reconstructed (by unfolding), the aspect ratio of both semi-axes is given by
#' the shape factor \code{S} between \eqn{(0,1]} and the orientation angle is named \code{alpha}.
#'
#' @param size	  matrix of axes lengths
#' @param angle   orientation angle of spheroids, see details
#' @param type    name of spheroid type, either "\code{prolate}" or "\code{oblate}"
#'
#' @return 		  section profiles object, either of class "\code{prolate}" or "\code{oblate}"
#' 
#' @author M. Baaske
#' @rdname sectionProfiles
#' @export
sectionProfiles <- function(size,angle,type=c("prolate","oblate")) {
	type <- match.arg(type)
	stopifnot(is.matrix(size))
	if(anyNA(size) || any(size<0))
		stop("'size' must have non-negative values.")
	if(anyNA(angle) || !is.numeric(angle) || any(angle<0))
		stop(paste("'angle' must have values between zero and ",quote(pi/2),sep=""))
	if(max(angle)>pi/2)
	 angle <- sapply(angle,.getAngle)	
	
    structure(list("A"=if(type=="prolate") size[,2] else size[,1],
				   "S"=size[,2]/size[,1],
				   "alpha"=angle),
		   class=type)
}

#' Simulation of a Poisson germ-grain process
#'
#' Simulation of a Poisson germ-grain process with either spheres, spheroids or spherocylinders as grains
#'
#' The function simulates a Poisson germ-grain process given the simulation parameters in the argument \code{theta}
#' and a predefined bounding box where the positions of the germs are generated independently according to a Poisson process with mean
#' intensity parameter \code{lam}. The function either generates \code{type="prolate"} or \code{type="oblate"} spheroids, spheres or spherocylinders.
#' The argument \code{size} is a name of the distribution function of the size of the grains, i.e. the semi-major axis lengths in case of spheroids,
#' radii for spheres or the lengths of the main axis of cylinders. 
#' 
#' The following directional orientation distributions of the spheroid's major-axis or cylinder's main axis are available:
#' a uniform distribution (\code{runifdir}), isotropic random planar distribution (\code{rbetaiso}, see the reference below)
#' and a "von Mises-Fisher" (\code{rvMisesFisher}) distribution. The simulation box is a list containing the ranges of each box dimension
#' corresponding to the lower and upper limits in each direction. If the argument \code{box} contains only a single range, i.e. \code{list(c(0,1)}, the
#' same limits are assumed for the remaining dimensions. The argument \code{rjoint} names a joint distribution function which can be any function
#' provided by the user in order to generate the required random parameters for spheroids or cylinders, see below.
#' 
#' In addition, the function supports an exact simulation [2] of the grains. In case of spheroids and cylinders setting \code{size="rbinorm"}
#' declares a bivariate normal size-shape distribution for which the exact simulation is available as an option if \code{perfect=TRUE}. More specifically,
#' for a bivariate normal vector \eqn{[X,Y]} with correlation parameter \eqn{\rho}, the length of the semi-major axis of a spheroid is given by \eqn{a=exp(x)}
#' with a (logit-transformed) shape parameter as \eqn{s=1/(1+exp(-y))} and thus a scaled semi-minor axis length \eqn{c=a*s}. This modification leads to a lognormally
#' distributed length of the semi-major axis. In case of cylinders, the lognormally distributed length of a cylinder is \eqn{len=h+2*r} where
#' \code{h} is the height and \eqn{r=len/2*s} the radius. The direction \code{u} of the a spheroid or cylinder is determined by the
#' major axis independent of the size and shape. Also, the following univariate distributions of the lengths \code{a} and \code{len}, the shapes \code{s} are
#' available: `\code{rbeta}`, `\code{rgamma}`, `\code{rlnorm}` and `\code{runif}`. Use `\code{const}` for a simulation with a constant length or shape.
#' Only simulations by "\code{rbinorm}" can use the exact simulation if \code{perfect=TRUE}.
#'
#' For spheres any distribution of the radii can be specified as a name of a user-defined function in \code{size} as long as the formal named
#' function parameters match the actual named parameters exactly as defined in the parameters given by \code{theta}.
#' Besides this, all other distributions given above are also available. In particular, setting \code{size="rlnorm"} leads to lognormally distributed
#' radii of the spheres in which case the exact simulation is available as an option (\code{perfect=TRUE}). Use "\code{const}" for a constant
#' radius of simulated spheres. 
#'  
#' The argument \code{pl>=0} denotes both the print level of intermediate information and type of return value after simulations. 
#' If \code{pl=10} then shorter lists of spheroids, spheres or cylinders are returned to speed up computation. Note that, the current
#' implementation does not include routines for unfolding the joint 3D size-shape-orientation
#' distribution of cylinders so far. 
#' 
#'
#' @param theta 	list of simulation parameters which must consist of elements: \code{size}, \code{shape} and \code{orientation}
#' @param lam   	mean number of grains per unit volume
#' @param size  	name of size distribution function 
#' @param shape 	name of shape distribution function  
#' @param orientation name of direction distribution function
#' @param type     spheroid type, either "\code{prolate}" or "\code{oblate}"
#' @param rjoint   name of joint distribution function given by the user
#' @param box	   simulation box
#' @param mu	   main orientation axis, \code{mu=c(0,0,1)} (default), as the alignment of the coordinate system
#' @param dz	   distance of the intersecting x-, y-plane to the origin
#' @param n		   normal vector of intersting plane
#' @param intersect use "\code{full}" to return the simulated system together with section profiles as two lists, choose "\code{only}" for section
#' 					 profiles only and "\code{original}" for the 3D system only
#' @param intern   logical, \code{intern=FALSE} (default), whether to return only section profiles with centers inside the simulation window
#' @param perfect  logical, \code{perfect=TRUE} (default) simulate perfect
#' @param pl  	   integer, print level and return value type, see details
#' @param label    character, label set to each generated grain, ´\code{N}´ (default)
#' 
#' @return 		   list of grains either of class "\code{prolate}" or "\code{oblate}", "\code{sphere}" or "\code{cylinder}" or a list of length two if
#' 				   \code{intersect="full"} where the first element stores the 3D generated objects as a list and the second the corresponding section profiles
#' 				   either in a short verion (if \code{pl=10}) or in a long version otherwise
#'
#' @example inst/examples/sim.R
#'
#' @references
#'	 \itemize{
#'		\item{} {Ohser, J. and Schladitz, K. 3D images of materials structures Wiley-VCH, 2009}
#'      \item{} {C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms.
#' 					Springer, Berlin, 2002. Zbl 0990.86007}
#' 	  }
#' @author M. Baaske
#' @rdname simPoissonSystem
#' @export
simPoissonSystem <- function(theta, lam, size="const", shape="const", orientation="rbetaiso",
								type=c("prolate","oblate","spheres","cylinders"), rjoint=NULL, box=list(c(0,1)),
								 mu=c(0,0,1), dz=0, n=c(0,1,0), intersect=c("original","full","only"), 
								  intern=FALSE, perfect=FALSE, pl=0, label="N")
{
	it <- pmatch(type,c("prolate","oblate","spheres","cylinders"))
	if(length(it)==0 || is.na(it))
	  stop(paste("`type` is one of: ", paste0("`",c("prolate","oblate","sphere","cylinder"),"`",collapse=","),sep=""))
	
	if(!is.numeric(lam) || !(lam>0) )
		stop("Expected 'lam' as non-negative numeric argument")

	if(length(box)==0 || !is.list(box))
	  stop("Expected argument 'box' as list type.")
	if(length(box)==1)
	  box <- rep(box[1],3)
  	else if(length(box)!=3)
	  stop("Simulation box has wrong dimensions.")
    names(box) <- c("xrange","yrange","zrange")

	# spheroid type
	type <- match.arg(type)
	intersect <- match.arg(intersect)
	
	cond <-
	 if(!is.null(rjoint))
	 {
		if(!is.function(rjoint))
		 stop("Unknown multivariate distribution function.")
		#largs <- theta[-(which(it==1))]
		if(!is.list(theta) || length(theta) == 0L)
		 stop("Expected 'theta' as list of named arguments matching exactly with formal arguments of function `rjoint`.")
				
		it <- match(names(theta),names(formals(rjoint)))
		funname <- as.character(substitute(rjoint))
		if(length(it) == 0L || anyNA(it))
			stop(paste("Arguments must match formal arguments of function `",funname,"`.",sep=""))

		# check function
		funret <- try(do.call(rjoint,theta))		
		if(inherits(funret,"try-error"))
		  stop(paste("Error in user defined function `",funname,"`.",sep=""))
		
	    fargs <- 
		  if(type %in% c("prolate","oblate")) {			  
			  if(!is.list(funret))
		  	     stop("Expected a list as return type in user defined function.")
			  c("a","b","c","u","shape","theta","phi")
		  } else if(type %in% c("cylinders")) {
			  if(!is.list(funret))
			     stop("Expected a list as return type in user defined function.")
			  c("h","r","u","theta","phi")		  
		  } else if(type %in% c("spheres")) {
			  if(!is.numeric(funret))
				  stop("Expected a numeric vector as return type in user defined function.")
			  c("r") 
		  } else { stop("Undefined type!")}
  
	    if(any(!(fargs %in% names(funret)))){		    
		 stop("Argument names of return value list does not match required arguments.")
	    }

		list("type"=type,"rdist"=funname,"box"=box,
			 "lam"=lam,"pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label,
			 "dz"=dz, "nsect"=n, "intern"=as.integer(intern),
			 "perfect"=as.integer(perfect), "intersect"=intersect)
			
	} else  {
		
		# set defaults if not provided
		if(!is.list(theta))
		 stop("Expected 'theta' as list of named arguments `size`, `shape`, `orientation`.")
		nms <- c("size","shape","orientation")
		it <- match(names(theta), nms)
		if(length(it)>0L && anyNA(it))
		 names(theta)[it] <- nms[it]		
				
		if(shape == "const")
		{
			if(!is.list(theta$shape) || length(theta$shape) == 0L)
			 theta$shape <- list(1,1)
		
		} else if(exists(shape, mode="function")) {
			# check arguments of supported distribution functions
			fargs <- names(formals(shape))
			if(shape %in% c("rbeta","rgamma","runif"))
				fargs <- fargs[-1]	
			if(!is.list(theta$shape) || length(theta$shape) == 0L)
			 stop("`shape` parameters must be given in `theta`.")
			it <- match(names(theta$shape),fargs)
			if(length(it)==0 || anyNA(it))
				stop(paste("Arguments of 'shape' must match formal arguments of function ",shape,sep=""))
		} else {
			stop(paste("Undefined `", shape, "` distribution function."))
		}
		
		if(!is.list(theta$size) || length(theta$size)==0L)
			stop("Arguments for size distribution must be given in `theta`.")
		if(size == "rbinorm") {
			 fargs <- c("mx","my","sdx","sdy","rho","kappa")
		     it <- match(names(theta$size),fargs)
			 if(length(it)==0 || anyNA(it))
				 stop(paste("Arguments of 'size' must match formal arguments: ", paste0("`",fargs,"`",collapse=",")))
			 
  	    } else if(exists(size, mode="function")) {
			 # check arguments of supported distribution functions
			 fargs <- names(formals(size))
			 if(size %in% c("rlnorm","rbeta","rgamma","runif"))
				 fargs <- fargs[-1]
			 
			 it <- match(names(theta$size),fargs)
			 if(length(it)==0 || anyNA(it))
				 stop(paste("Arguments of 'size' must match formal arguments of function ",size,sep=""))
			 		
		 } else {
			stop(paste("Undefined `", size, "` distribution function."))
		 }
		 
		 it <- pmatch(orientation,c("runifdir","rbetaiso","rvMisesFisher"))
		 if(is.na(it) && !exists(orientation, mode="function"))
			 stop("Undefined distribution function for the orientation/direction.")	 
		 if(!is.list(theta$orientation) || length(theta$orientation) == 0L)
			 theta$orientation <- list("kappa"=1)
		 
		 list("type"=type, "lam"=lam,
			  "rdist"=list("size"=size, "shape"=shape, "orientation"=orientation),
			  "box"=box,"pl"=pl,"mu"=mu,"rho"=.GlobalEnv,
			  "dz"=dz, "nsect"=n, "intern"=as.integer(intern),"label"=label,
			  "perfect"=as.integer(perfect), "intersect"=intersect)
	}
		 
	.Call(C_PoissonSystem, theta, cond)	
}

#' Calculate coefficients for unfolding
#'
#' Calculate coefficients of discretized integral equation
#'
#' In order to apply the EM algorithm to the stereological
#' unfolding procedure for the joint size-shape-orientation distribution
#' one first has to calculate the coefficients of the discretized integral
#' equation. This step is the most time consuming part of unfolding
#' the parameters and therefore has been separated in its own function.
#' Further, the number of classes for size, shape and orientation do not
#' need to be equal, whereas the same class limits are used for binning spatial and planar values.
#' This might be changed in future releases.
#' Using multiple cpu cores is controlled by either setting the option \code{par.unfoldr} in the global
#' R environment or passing the number of cores \code{nCores} directly.
#'
#' @param   breaks  list of bin vectors
#' @param   stype   either \code{prolate} or \code{oblate}
#' @param   check   logical, whether to run some input checks
#' @param   nCores  number of cores used to calculate the coefficients
#' @return  coefficient array
#'
#' @example inst/examples/coeffarray.R
#' 
#' @author M. Baaske
#' @rdname coefficientMatrixSpheroids 
#' @export
coefficientMatrixSpheroids <- function(breaks, stype=c("prolate","oblate"),
										check=TRUE,nCores=getOption("par.unfoldr",1))
{
	stype <- match.arg(stype)
	it <- match(names(breaks), c("size","angle","shape"))
	if (length(it)==0 || anyNA(it))
		stop("Expected 'breaks' as named list of: 'size','angle','shape' ")
	if (is.unsorted(breaks$size) || is.unsorted(breaks$angle) || is.unsorted(breaks$shape))
		stop("'breaks' list must contain non-decreasingly sorted classes")

	if(check) {
		if(any(breaks$size<0))
			stop("Breaks vector 'size' must have non-negative values.")
		if(min(breaks$angle)<0 || max(breaks$angle)>pi/2)
			stop(paste("Breaks vector 'angle' must have values between zero and ",pi/2,sep=""))
		if(min(breaks$shape)<0 || max(breaks$shape)>1)
			stop("Breaks vector 'shape' must have values between 0 and 1.")
	}

	.Call(C_CoefficientMatrixSpheroids,
			breaks$size,breaks$angle,breaks$shape,
			breaks$size,breaks$angle,breaks$shape, list(stype,nCores))

}

#' Spheroid vertical section
#'
#' Get a vertical section of a spheroid system
#'
#' The function intersects a spheroid system by an intersecting plane defined by a normal 
#' vector \code{n=c(0,1,0)} (default), which depends on the main orientation axis of the
#' coordinate system.
#'
#' @param S		 list of spheroids, see \code{\link{simPoissonSystem}}
#' @param d 	 distance of the box-aligned intersecting plane from the origin
#' @param n 	 normal vector which defines the intertecting plane
#' @param intern logical, \code{FALSE} (default), return all section profiles otherwise
#' 				 only those which have their centers inside the intersection window (if the 
#' 				 intersected spheroid system had been simulated using exact simulation)
#' @return 	 	 list sizes \code{A}, shapes \code{S} and angles \code{alpha} of section profiles
#' 
#' @author M. Baaske
#' @rdname verticalSection 
#' @export
verticalSection <- function(S,d,n=c(0,1,0),intern=FALSE) {
	if(!(class(S) %in% c("prolate","oblate")))
	 stop("Only spheroids of class `prolate` or `oblate` can be used.")
 	stopifnot(is.numeric(d))	
 	stopifnot(is.logical(intern))
	if( sum(n)>1 )
	  stop("Normal vector is like c(0,1,0). ")
	
    ss <- Call(C_IntersectPoissonSystem,
		 		 as.character(substitute(S)),
		  		 list("nsect"=n,"dz"=d,"intern"=intern,"pl"=10),
		  		 .GlobalEnv)
  	
	A <- if(class(S)=="prolate")
			sapply(ss,"[[",2)
		 else sapply(ss,"[[",1)
	
    structure(
	    list("A"=A,
			 "S"=sapply(ss,"[[",3),
		 "alpha"=sapply(ss,"[[",4)),
	  class=class(S)
	)
}

#' Intersect a 3D system of grains
#'
#' Intersect a system of spheres, spheroids or cylinders by a given plane
#'
#' The function intersects a given system of spheres, spheroids or cylinders by a plane which is defined by a
#' normal vector, e.g. \code{n=c(0,1,0)}, perpendicular to one of the lateral planes of the simulation box. 
#' The print level \code{pl>=0} sets the type of return value. If \code{pl=10} the functions only returns
#' the lengths of the semi axes, the shape factor of the angle in the intersection plane between \code{[0,pi/2]}.    
#'
#' @param S		 list of spheroids, see \code{\link{simPoissonSystem}}
#' @param d 	 distance of the the box-aligned intersecting plane from the origin
#' @param n 	 normal vector which defines the intersting plane
#' @param intern logical, \code{FALSE} (default), return all section profiles otherwise
#' 				 only those which have their centers inside the intersection window (if the 
#' 				 intersected spheroid system had been simulated using exact simulation)
#' @param pl	 integer, \code{pl=0} (default), only return pointer to stored intersections
#' 				
#' @return list of size, shape and angle of section profiles
#' 
#' @author M. Baaske
#' @rdname intersectSystem
#' @export
intersectSystem <- function(S, d, n=c(0,1,0), intern=FALSE, pl=0) {
	if(!(class(S) %in% c("prolate","oblate","spheres","cylinders")))
	  stop("Unknown object class.")
    stopifnot(is.numeric(d))
	stopifnot(is.logical(intern))
	if(sum(n)>1)
	  stop("Normal vector is like c(0,1,0). ")	
	
	.Call(C_IntersectPoissonSystem,
			as.character(substitute(S)),
			list("nsect"=n,"dz"=d,"intern"=intern,"pl"=pl),
		    .GlobalEnv)	
}

#' Generate 3D plot of spheroid system
#'
#' Draw a spheroid system
#'
#' The function requires the package \code{rgl} to be installed.
#'
#' @param S				a list of spheroids
#' @param box			simulation box
#' @param draw.axes		logical: if true, draw the axes
#' @param draw.box	    logical: if true, draw the bounding box
#' @param draw.bg	    logical: if true, draw the a gray background box
#' @param bg.col 		background color used to draw the box background
#' @param clipping 		logical: if true clip to the bounding box
#' @param ...			further material properties passed to 3d plotting functions
#' 
#' @return NULL
#' @author M. Baaske
#' @rdname spheroids3d
#' @export
spheroids3d <- function(S, box, draw.axes=FALSE, draw.box=TRUE, draw.bg=TRUE, bg.col="white", clipping=FALSE, ...)
{
	if (!requireNamespace("rgl", quietly=TRUE))
	 stop("Please install 'rgl' package from CRAN repositories before running this function.")

	ellipsoid3d <- function(rx=1,ry=1,rz=1,n=50,ctr=c(0,0,0), qmesh=FALSE,trans = rgl::par3d("userMatrix"),...) {
		if (missing(trans) && !rgl::rgl.cur())
			trans <- diag(4)
		degvec <- seq(0,2*pi,length=n)
		ecoord2 <- function(p) {
			c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2]))
		}
		v <- apply(expand.grid(degvec,degvec),1,ecoord2)
		if (qmesh)
			v <- rbind(v,rep(1,ncol(v))) ## homogeneous
		e <- expand.grid(1:(n-1),1:n)
		i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
		i2 <- i1+1
		i3 <- (i1+n-1) %% n^2 + 1
		i4 <- (i2+n-1) %% n^2 + 1
		i <- rbind(i1,i2,i4,i3)
		if (!qmesh)
			rgl::quads3d(v[1,i],v[2,i],v[3,i],...)
		else
			return(rgl::rotate3d(rgl::qmesh3d(v,i,material=...),matrix=trans))
	}
	sphere <- ellipsoid3d(qmesh=TRUE,trans=diag(4))

	spheroid3d <- function (x=0,y=0,z=0, a=1, b=1,c=1, rotM, subdivide = 3, smooth = TRUE){
		result <- rgl::scale3d(sphere, a,b,c)
		result <- rgl::rotate3d(result,matrix=rotM)
		result <- rgl::translate3d(result, x,y,z)
		invisible(result)
	}
	N <- length(S)
	ll <- lapply(S, function(x)
				spheroid3d(x$center[1],x$center[2],x$center[3],
						x$acb[1],x$acb[2],x$acb[3],rotM=x$rotM))

	rgl::shapelist3d(ll,...)

	x <- box$xrange[2]
	y <- box$yrange[2]
	z <- box$zrange[2]
	## draw gray box
	if(draw.bg) {
		c3d.origin <- rgl::translate3d(rgl::scale3d(rgl::cube3d(col=bg.col, alpha=0.1),x/2,y/2,z/2),x/2,y/2,z/2)
		rgl::shade3d(c3d.origin)
	}
	if(clipping) {
		rgl::clipplanes3d(-1,0,0,x)
		rgl::clipplanes3d(0,-1,0,y)
		rgl::clipplanes3d(0,0,-1,z)
		rgl::clipplanes3d(1,0,0,0)
		rgl::clipplanes3d(0,1,0,0)
		rgl::clipplanes3d(0,0,1,0)

	}

	if(draw.axes) {
		rgl::axes3d(c('x','y','z'), pos=c(0,0,0))
		rgl::title3d('','','x','y','z')
	}
	## draw box
	if(draw.box) {
		rgl::axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
				expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
	}
}


#' 3D a cylinder system
#'
#' Draw 3D spherocylinders
#'
#' The function requires the package \code{rgl} to be installed.
#'
#' @param S				a list of cylinders
#' @param box			simulation box
#' @param draw.axes		logical, if \code{TRUE}, draw the axes
#' @param draw.box	    logical, if \code{TRUE}, draw the bounding box
#' @param clipping 		logical, if \code{TRUE}, clip to the bounding box
#' @param ...			further material properties passed to 3d plotting functions
#' 
#' @return NULL
#' @author M. Baaske
#' @rdname cylinders3d
#' @export
cylinders3d <- function(S, box, draw.axes=FALSE, draw.box=TRUE, clipping=FALSE,...)
{
	if (!requireNamespace("rgl", quietly=TRUE))
		stop("Please install 'rgl' package from CRAN repositories before running this function.")
	
	cylinder <- function(m, radius=1, h=1, rotM=diag(3), u=c(0,0,1)) {
		cyl <- rgl::cylinder3d(rbind(c(0,0,0),c(0,0,1)), radius=radius,
				e1=cbind(0, 0, 1), e2=cbind(1, 0, 0), sides=25,closed=-2	)
		result <-rgl::scale3d(cyl,1,1,h)
		result <- rgl::rotate3d(result,matrix=rotM)
		## translate to midpoint of cylinder, center of mass
		m <- m-0.5*h*u
		result <- rgl::translate3d(result, m[1],m[2],m[3])
		invisible(result)
	}
	
	args <- list(...)
	ok <- sapply(S, function(x) x$length>0)
	
	cyls <- lapply(S[ok], function(x) { cylinder(x$center, x$r, x$length, x$rotM, x$u) })
	rgl::shapelist3d(cyls,...)
	
	if("col" %in% names(args)) {
		cols <- rep(rep(args$col,length.out=length(S)),each=2)
		args$col <- NULL
	} else cols <- "black"
	
	# spheres
	Xc <- do.call(rbind,lapply(S[ok],function(x) rbind( c(x$origin0,x$r),c(x$origin1,x$r))))
	rgl::spheres3d(Xc,radius=Xc[,4],col=cols,unlist(args))
	if(!all(ok)) {
		Xc <- do.call(rbind,lapply(S[!ok],function(x) rbind(c(x$center,x$r))))
		rgl::spheres3d(Xc,radius=Xc[,4],col="darkgray",unlist(args))
	}
	
	x <- box$xrange[2];	y <- box$yrange[2];	z <- box$zrange[2]
	c3d.origin <- rgl::translate3d(rgl::scale3d(rgl::cube3d(col="darkgray", alpha=0.1),x/2,y/2,z/2),x/2,y/2,z/2)
	rgl::shade3d(c3d.origin)
	
	if(clipping) {
		rgl::clipplanes3d(-1,0,0,box$xrange[1])
		rgl::clipplanes3d(0,-1,0,box$yrange[1])
		rgl::clipplanes3d(0,0,-1,box$zrange[1])
		rgl::clipplanes3d(1,0,0,box$xrange[2])
		rgl::clipplanes3d(0,1,0,box$yrange[2])
		rgl::clipplanes3d(0,0,1,box$zrange[2])
	}
	
	if(draw.axes) {
		rgl::axes3d(c('x','y','z'), pos=c(0,0,0))
		rgl::title3d('','','x','y','z')
	}
	if(draw.box) {
		rgl::axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
				expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
	}
}


#' Plot Spheroid intersection
#'
#' Drawing section profiles in 3D plane
#'
#' The function requires the package \code{rgl} to be installed.
#'
#' @param E				a list of spheroid intersections
#' @param n 			the normal vector of the intersecting plane
#' @param np			number of points for polygon approximation of ellipses
#' @return NULL
#' @author M. Baaske
#' @rdname drawSpheroidIntersection
#' @export
drawSpheroidIntersection <- function(E, n=c(0,1,0), np=25) {
	ind <- which(n==0)
	.pointsOnEllipse <- function(E,t) {
		E$phi <- E$phi + 0.5*pi
		xt <- E$center[1] + E$ab[1]*cos(t)*cos(E$phi)-E$ab[2]*sin(t)*sin(E$phi)
		zt <- E$center[2] + E$ab[1]*cos(t)*sin(E$phi)+E$ab[2]*sin(t)*cos(E$phi)
		yt <- rep(0,length(t))
		m <- matrix(0,nrow=length(xt),ncol=3)
		m[,ind[1]] <- xt
		m[,ind[2]] <- zt
		m
	}
	.plotEllipse <- function(x) {
		M <- .pointsOnEllipse(x,t)
		rgl::polygon3d(M[,1],M[,2],M[,3],fill=TRUE,coords=ind)
	}
	s <- 2*pi/np
	t <- seq(from=0,to=2*pi,by=s)
	invisible(lapply(E,function(x) .plotEllipse(x) ))
}



.checkArgs <- function(optlist, options) {
	if(is.null(names(options)))
		stop("Options should be a list of named arguments.")
	if (!is.list(options) || "" %in% names(options))
		stop("Argument ",as.character(substitute(options)), "  must be a list of named (character) elents.")
	optnames <- (names(options) %in% names(optlist))
	if (!all(optnames)) {
		unames <- as.list(names(options)[!(optnames)])
		stop(paste(c("Unknown arguments in ",as.character(substitute(options))," : ",do.call("paste", c(unames, sep = ", "))), collapse=" "))
	}	
	return (0)
}

#' Sphere planar section
#'
#' Intersection of given spheres by a plane
#'
#' The function intersects a sphere system and returns the diameters of the corrsponding planar section profiles.
#'
#' @param S 		list of spheres of class \code{sphere}, see \code{\link{simPoissonSystem}}
#' @param d 		distance of the intersecting xy-plane to the origin (as a planar section)
#' @param intern 	logical, \code{FALSE} (default), return all planar section profiles otherwise
#' 					only those which have their centers inside the intersecting window
#' @param pl		print level, \code{pl=0} (default) for no output
#'
#' @return 			numeric vector of planar section diameters
#' 
#' @author M. Baaske
#' @rdname planarSection
#' @export
planarSection <- function(S, d, intern=FALSE, pl=0) {
	stopifnot(is.logical(intern))
	if(!(c("spheres") %in% class(S) ))
		stop("Expected spheres as list argument.")
	
	sp <- .Call(C_IntersectPoissonSystem,
			as.character(substitute(S)),
			list("nsect"=c(0,0,1),"dz"=d,"intern"=intern,"pl"=pl),
			.GlobalEnv)
	
	if(is.list(sp))								# full list of section profiles
		return (sapply(sp,function(x) 2.0*x$r))		
	else return (2*sp)							# only diameters are returned
}

#' Binning numeric values
#'
#' Vector of count data
#'
#' This function provides basic binning (grouping) of numeric values into
#' classes defined by the breaks vector \code{bin}. The values are binned
#' according to bin[i[j]]\eqn{<}x[j]\eqn{\leq} bin[i[j]+1] for interval i=1,...,N-1
#' for \code{length(bin)=N} and value x[j].
#' If x[j] > bin[N] or x[j] < bin[1] then x[j] is not counted at all.
#'
#' @param x   	 numeric values to be binned
#' @param bin    non-decreasingly sorted breaks vector
#' @param na.rm  logical, removing missing values (including NaN) in the argument \code{x}?
#'
#' @return Vector of count data
#'
#' @examples
#' 	x <- runif(100,0,1)
#' 	bin <- seq(0,1,by=0.1)
#' 	binning1d(x,bin)
#' 
#' @author M. Baaske
#' @rdname binning1d
#' @export
binning1d <- function(x,bin, na.rm = FALSE) {
	if (anyNA(x)) {
		if(na.rm) x <- x[which (!is.na(x))]
		else stop("Vectors contains missing values or NAs.")
	}
	if(anyNA(bin) )
		stop("'bin' vector contains missing values or NAs.")
	
	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")
	
	.Call(C_Binning1d,x,bin)
}

#' Calculate coefficients (spheres)
#'
#' Matrix of coefficients for Wicksell's corpuscle problem
#'
#' The function calculates the matrix of coefficients of the
#' discretized integral equation for Wicksell's corpuscle problem.
#'
#' @param bin  non-decreasing vector of class limits
#' @return 	   array of coefficients
#'
#' @references
#'  Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
#' 
#' @author M. Baaske
#' @rdname coefficientMatrixSpheres
#' @export
coefficientMatrixSpheres <- function(bin) {
	## count of bin classes
	n <- length(bin)-1
	if (n<=0)
		stop("'bin' has length zero")
	if (anyNA(bin))
		stop("Vectors contain NA values")
	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")
	
	if(any(bin<0))
		stop("Breaks vector must have non-negative values.")
	
	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(bin),
			p=as.numeric(matrix(0,n,n)))$p
	dim(p) <- c(n,n)
	return (p)
}

#' Expectation Maximization algorithm
#'
#' Estimation of empirical sphere diameter distribution
#'
#' The function performs the EM algorithm.
#'
#' @param y   		vector of observed absolute frequencies of circle diameters
#' @param bin 		non-decreasing vector of class limits
#' @param maxIt		maximum number of iterations used
#' @return 			vector of count data of absolute frequenties of sphere diameters
#'
#' @example inst/examples/sphere.R
#'
#' @references
#' Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
#' 
#' @author M. Baaske
#' @rdname em.saltykov
#' @export 
em.saltykov <- function(y,bin,maxIt=32) {
	if (length(y)==0) stop("input array 'y' has length zero")
	if (anyNA(y) || anyNA(bin) )
		stop("Vectors contain NA values")
	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")
	
	n <- length(y)
	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(bin),
			p=as.vector(matrix(0,n,n)))$p
	
	theta <- y+1.0e-6
	.C(C_em_saltykov,as.integer(n), as.integer(maxIt),
			as.numeric(p),as.numeric(y),theta=as.numeric(theta))$theta
}

#' Digitize section profiles
#'
#' Digitize 2D objects (section profiles)  
#'
#' A list of section profiles, e.g. discs, ellipses or segments from intersected spherocylinders, is digitized according
#' to a given resolution such that the result is an image stored in matrix form.  
#'
#' @param sp   		list of section profiles, see \code{\link{intersectSystem}}
#' 
#' @return 			image as a matrix
#'
#' @author M. Baaske
#' @rdname digitizeProfiles
#' @export 
digitizeProfiles <- function(sp, delta, win = NULL ) {
	stopifnot(is.numeric(delta))
	if(!is.null(win)){
		if(!is.list(win) || length(win) == 0L )
			stop("Expected argument 'win' as list type.")
		if(length(win) == 1L)
			win <- rep(win[1],2)
		else if(length(win) != 2L)
			stop("Digitization window has wrong dimensions.")		
	}
	if(delta > 1.0)
	  stop("´delta` should be equal or less than 1.")
  
	.Call(C_DigitizeProfiles,as.character(substitute(sp)), delta, win, .GlobalEnv)
}