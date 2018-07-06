###############################################################################
# Author:  M. Baaske
# Date:	   2018/06/15	
# File:    cylinder.R: 
# 
# Comment: Simulation of sphero cylinders, intersection and visualization
# 
###############################################################################

#' Simulation of cylinder system
#'
#' Simulation of Poisson cylinder system
#'
#' The function simulates a Poisson (sphero)cylinder system according to the supplied
#' simulation parameter \code{theta} in a predefined simulation box.
#' The argument \code{size} is of type string and denotes the major-axis length random generating
#' function name.
#'
#' For the directional orientation of the cylinder axis one has the choice of a uniform
#' (\code{runifdir}), isotropic random planar (\code{rbetaiso}, see reference) or von Mises-Fisher
#' (\code{rvMisesFisher}) distribution. The simulation box is a list containing of vector arguments
#' which correspond to the lower and upper points in each direction. If the argument \code{box} has
#' only one element, i.e. \code{list(c(0,1)}, the same extent is used for the other dimensions.
#' If \code{rjoint="rmulti"} names a joint random generating function then argument \code{size} is ignored
#' (see example file "sim.R").
#' For the purpose of exact simulation [2] setting \code{size="rbinorm"} declares a bivariate normal
#' size-shape distribution. For a bivariate normal vector \eqn{[X,Y]} with correlation parameter \eqn{\rho}
#' the length of the cylinder is defined as \eqn{exp(x)} with (logit-transformed) shape parameter \eqn{s=1/(1+exp(-y))}.
#' This modification leads to a lognormally distributed length \code{h+2*r} of the cylinder where \code{h} is the height 
#' and \code{r} the radius of the (sphero)cylinder (and also the caps). The direction axis \code{u} is along the longer side of the
#' cylinder and independent of the size and shape. The following univariate distributions of the cylinder length and shape are also available:
#' `\code{rbeta}`, `\code{rgamma}`, `\code{rlnorm}` and `\code{runif}` distribution. Use `\code{const}` for a constant length or shape (radius) of the cylinders.
#' Only for `\code{rbinorm}` the exact simulation type can be used as an option if \code{perfect=TRUE}.
#' The current implementation does not include routines for unfolding the joint 3D size-shape-orientation
#' distribution of cylinders so far.
#'
#' The argument \code{pl} denotes the print level of output information during simulation.
#' Currently, only \code{pl}=0 for no output and \code{pl}>100 for some additional info is implemented.
#'
#' @param theta 		simulation parameters (see examples)
#' @param lam   		mean number of spheroids per unit volume
#' @param size  		\code{size="const"} (default) or name of random generating function
#' 						a for specific size distribution
#' @param shape 		either \code{shape="const"} as a constant portion of the axis length
#' 					    or \code{shape="rbeta"} as beta distributed shape factor. For a radius distribution
#'	 					use your own joint distribution, see details.
#' @param orientation   name of random generating function for orientation distribution
#' @param rjoint 		name of joint random generating function
#' @param box 			simulation box
#' @param mu  			main orientation axis, \code{mu=c(0,1,0)} (default)
#' @param dz		    distance of the intersecting x-, y-plane to the origin
#' @param n			    normal vector of intersting plane
#' @param profiles 		logical, \code{profiles=FALSE} (default), whether the simulated system of cylinder is intersected afterwards in which case
#' 					    only sections profiles are returned 
#' @param intern        logical, \code{intern=FALSE} (default), whether to return only section profiles with centers inside the simulation window
#' @param perfect 		logical, \code{perfect=TRUE} (default) simulate perfect
#' @param pl  			integer, print level and return value definition
#' @param label 		character, label passed to each simulated cylinder, `\code{N}` (default)
#' @return 				list of cylinders
#'
#' @example inst/examples/cylinder.R
#'
#' @references
#'	 \itemize{
#'		\item{} {Ohser, J. and Schladitz, K. 3D images of materials structures Wiley-VCH, 2009}
#'      \item{} {C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms.
#' 					Springer, Berlin, 2002. Zbl 0990.86007}
#' 	  }
#' @author M. Baaske
#' @rdname simCylinderSystem
#' @export
simCylinderSystem <- function(theta, lam, size="const", shape="const",
	orientation="rbetaiso", rjoint=NULL, box=list(c(0,1)), mu=c(0,1,0),
	dz=0, n=c(0,1,0), profiles = FALSE, intern=FALSE, perfect=TRUE, pl=0, label="N")
{
	if(!is.list(theta))
		stop("Expected 'theta' as list of named  arguments.")
	if(any(sapply(theta,length) == 0L))
		stop("At least one of the lists of simulation parameters given in 'theta' has length zero.")
	if(!is.numeric(lam) || !(lam>0) )
		stop("Expected 'lam' as non-negative numeric argument")

	if(length(box)==0 || !is.list(box))
		stop("Expected argument 'box' as list type.")
	if(length(box)==1)
		box <- rep(box[1],3)
	else if(length(box)!=3)
		stop("Simulation box has wrong dimensions.")
	names(box) <- c("xrange","yrange","zrange")

	if(!is.null(rjoint)) {
		if(!exists(rjoint, mode="function"))
			stop("Unknown multivarirate random generating function.")
		#largs <- theta[-(which(it==1))]
		it <- match(names(theta),names(formals(rjoint)))
		if(length(it)==0 || anyNA(it))
			stop(paste("Arguments must match formal arguments of function ",rjoint,sep=""))

		# check function
		funret <- try(do.call(rjoint,theta))
		if(!is.list(funret))
			stop("Expected list as return type in user defined function.")
		if(inherits(funret,"try-error"))
			stop(paste("Error in user defined function ",rjoint,".",sep=""))
		if(any(!(c("h","r","u","theta","phi") %in% names(funret))))
			stop("Argument names of return value list does not match required arguments.")

		cond <- list("rdist"=rjoint,"box"=box,"perfect"=as.integer(perfect),
				"lam"=lam, "pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=as.character(label),
				"dz"=dz, "nsect"=n, "intern"=as.integer(intern))
		
	} else  {
	
		it <- match(names(theta), c("size","shape","orientation"))
		if(!is.list(theta) || anyNA(it))
			stop("Expected 'theta' as list of named arguments `size`, `shape`, `orientation`.")
		if(!is.list(theta$size) || !is.list(theta$shape) || !is.list(theta$orientation) )
			stop("Expected `size`, `shape`, `orientation` as lists.")
		
		it <- pmatch(orientation,c("runifdir","rbetaiso","rvMisesFisher"))
		if(is.na(it) && !exists(orientation, mode="function"))
			stop("Undefined distribution function for orientation/direction.")
		
		cond <- list("rdist"=list("size"=size,"shape"=shape,"orientation"=orientation),
					 "lam"=lam, "box"=box, "pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label,
					 "dz"=dz, "nsect"=n, "intern"=as.integer(intern),
					 "perfect"=as.integer(perfect))

		if(cond$rdist$shape != "const") {
			if(length(theta$shape)==0L)
				stop("Arguments for shape distribution must be given.")
		} else if(exists(cond$rdist$shape, mode="function")) {
				 # check arguments of supported distribution functions
				 fargs <- names(formals(cond$rdist$shape))
				 if(cond$rdist$shape %in% c("rbeta","rgamma","runif"))
					 fargs <- fargs[-1]				
				 it <- match(names(theta$shape),fargs)
				 if(length(it)==0 || anyNA(it))
					 stop(paste("Arguments of 'shape' must match formal arguments of function ",cond$rdist$shape,sep=""))
		} else
		    stop(paste("Undefined `", cond$rdist$shape, "` distribution function."))
					 
			 
		if(cond$rdist$size == "const") {
			 if(length(theta$size)==0L)
				 stop("Arguments for size distribution must be given.")
		} else if(cond$rdist$size == "rbinorm") {
			 fargs <- c("mx","my","sdx","sdy","rho","kappa")
			 it <- match(names(theta$size),fargs)
			 if(length(it)==0 || anyNA(it))
				 stop(paste("Arguments of 'size' must match formal arguments: ", paste0("`",fargs,"`",collapse=",")))
			 
		} else if(exists(cond$rdist$size, mode="function")) {
			fargs <- names(formals(cond$rdist$size))
			if(cond$rdist$size %in% c("rlnorm","rbeta","rgamma","runif"))
				fargs <- fargs[-1]
			
			it <- match(names(theta$size),fargs)
			if(length(it)==0 || anyNA(it))
				stop(paste("Arguments of 'size' must match formal arguments of function `",cond$rdist$size, "`.",sep=""))
			
		} else {
			stop(paste("Undefined `", cond$rdist$size, "`  distribution function."))
		}		
	}
	
	if(profiles) {
		.Call(C_SimulateCylindersAndIntersect, theta, cond)	
	} else {
		structure(.Call(C_CylinderSystem, theta, cond),
				"mu"= mu, "lam"=lam, "box" = box, "perfect"=perfect)
	}
	
}

.cylVol <- function(X) {
	sum(sapply(X, function(x) { pi*x$r^2*x$length + 4/3*pi*x$r^3 }))
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
cylinders3d <- function(S, box, draw.axes=FALSE, draw.box=TRUE, clipping=FALSE,...) {
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

#' Spherocylinder intersection
#' 
#' Intersect a spherocylinder system 
#' 
#' The function intersects a cylinder system by an intersecting plane defined by a
#' normal vector \code{n=c(0,1,0)} (default), which depends on the main orientation
#' axis of the coordinate system.
#'
#' @param S		 list of cylinders, see \code{\link{simCylinderSystem}}
#' @param d 	 distance of intersecting plane to the origin
#' @param n 	 normal vector of intersting plane
#' @param intern \code{intern=FALSE} (default) return all section profiles otherwise
#' 				 only those which have their centers inside the intersection window
#' @param pl	 integer, \code{pl=0} (default), other options are not yet available
#' 
#' @return list of size, shape and angle of section profiles
#' @author M. Baaske
#' @rdname cylinderIntersection
#' @export
cylinderIntersection <- function(S, d, n = c(0,1,0), intern=FALSE, pl=0) {
	stopifnot(is.logical(intern))
	if(sum(n) > 1)
		stop("Normal vector is like c(0,1,0). ")
	if(!(class(S) %in% c("cylinder")))
		stop("Class must be `cylinder`.")
	.Call(C_IntersectCylinderSystem, as.character(substitute(S)), n, d, intern, .GlobalEnv, pl)	
}


#' Cylinder intersection
#' 
#' Simulate a cylinder system and intersect
#' 
#' The function first simulates a cylinder system according to the parameter \code{theta} and only
#' returns the resulting section profiles.
#' 
#' @param theta simulation parameters
#' @param cond  conditioning object for simulation and intersection
#' 
#' @return list of intersection profiles
#' @author M. Baaske
#' @rdname simCylinderIntersection
#' @export
simCylinderIntersection <- function(theta, cond) {
	.Call(C_SimulateCylindersAndIntersect,c("lam"=cond$lam,theta), cond)
}
