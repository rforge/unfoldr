###############################################################################
# Author:  M. Baaske
# Date:	   09.05.2016
# File:    cylinder.R:
#
# Comment:
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
#' For the purpose of exact simulation setting \code{size} equal to \code{rbinorm} declares a bivariate normal
#' size-shape distribution which leads to a lognormally distributed half height (length) \code{h/2} of the
#' cylinder. The main orientation axis of the cylinder is called \code{u} where its length equals \code{h}
#' without the end caps and a scaled radius \code{r}. The total length then equals \code{h+2r}. If \eqn{[X,Y]} follow a bivariate normal distribution with 
#' correlation parameter \eqn{\rho} then \eqn{h=2.0*exp(x)} defines the sample cylinder axis length together
#' with the scaled radius \eqn{r=0.5*h*s} and shape parameter set to \eqn{s=1/(1+exp(-y))}. The parameter 
#' \eqn{\rho} defines the degree of correlation between the cylinder axis length and cylinder radius which 
#' must be provided as part of the list of simulation parameters \code{theta}. The method of exact simulation
#' is tailored to the above described model. For a general approach please see the given reference below.
#' Other (univariate) cylinder axis lengths types include the beta, gamma, lognormal and uniform distribution
#' where the shape factor to get the radius either follows a beta distribution or is set to a constant.
#' Despite the case of constant size simulations all other simulations are done as perfect simulations.
#' The current implementation does not include routines for unfolding the joint 3d size-shape-orientation
#' distribution of cylinders so far. However, this feature this might be provided in a later version.
#'
#' The argument \code{pl} denotes the print level of output information during simulation.
#' Currently, only \code{pl}=0 for no output and \code{pl}>100 for some additional info is implemented.
#'
#' @param theta 		simulation parameters
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
#' @param perfect 		logical: \code{perfect=TRUE} (default) simulate perfect
#' @param pl  			optional: print level
#' @param label 		optional: set a label to all generated spheroids
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
simCylinderSystem <- function(theta, lam, size="const", shape="const",
						orientation="rbetaiso", type=c("sphero","elong"), 
						 rjoint=NULL, box=list(c(0,1)), mu=c(0,1,0), perfect=TRUE, pl=0, label="N")
{
	it <- pmatch(type,c("sphero","elong"))
	if(length(it)==0 || is.na(it)) stop("Cylinder type 'type' must be either 'sphero' or 'elong'.")
	
	if(!is.list(theta))
		stop("Expected 'theta' as list of named  arguments.")
	if(!is.numeric(lam) || !(lam>0) )
		stop("Expected 'lam' as non-negative numeric argument")

	if(length(box)==0 || !is.list(box))
		stop("Expected argument 'box' as list type.")
	if(length(box)==1)
		box <- rep(box[1],3)
	else if(length(box)!=3)
		stop("Simulation box has wrong dimensions.")
	names(box) <- c("xrange","yrange","zrange")

	type <- match.arg(type)
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
		if(any(!(c("a","b","u","shape","theta","phi") %in% names(funret))))
			stop("Argument names of return value list does not match required arguments.")

		structure(.Call(C_CylinderSystem,
						list("lam"=lam,"rmulti"=theta),
						list("type"=type,"rdist"=rjoint,"box"=box,"perfect"=0,
							 "pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label)),
				box = box)

	} else  {
		theta <- c("lam"=lam,theta)
		it <- match(names(theta), c("lam","size","shape","orientation"))
		if(!is.list(theta) || anyNA(it))
			stop("Expected 'theta' as list of named arguments.")
		if(!is.list(theta$size) || !is.list(theta$shape) || !is.list(theta$orientation) )
			stop("Expected 'size','shape' and 'orientation' as lists of named arguments.")
		it <- pmatch(orientation,c("runifdir","rbetaiso","rvMisesFisher"))
		if(is.na(it) && !exists(orientation, mode="function"))
			stop("Unknown random generating function for orientation distribution.")

		if (missing(size))
		 stop("Argument 'size' has to be given if 'rjoint' is 'NULL'!")
	 	sdistr <- c("const","rbeta","rgamma","runif")
	 	its <- pmatch(shape,sdistr)
	 	if(length(its)==0 || is.na(its))
		 stop("Unknown shape distribution set. Only 'const', 'rbeta' supported.")

		cond <- list("type"=type,"rdist"=list("size"=size, "shape"=shape,"orientation"=orientation),
					 "box"=box, "pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label,"perfect"=as.integer(perfect))

		if(cond$rdist$size %in% c("const","rbinorm")) {
			structure(.Call(C_CylinderSystem, theta, cond), box = box)
		} else if(exists(cond$rdist$size, mode="function")) {
			fargs <- names(formals(cond$rdist$size))
			if(cond$rdist$size %in% c("rlnorm","rbeta","rgamma","runif"))
				fargs <- fargs[-1]

			it <- match(names(theta$size),fargs)
			if(length(it)==0 || anyNA(it))
				stop(paste("Arguments of 'size' must match formal arguments of function ",cond$rdist$size,sep=""))

			structure(.Call(C_CylinderSystem, theta, cond), box = box)
		} else
			stop(paste("The ", cond$rdist$size, "random generating function must be defined"))

	}
}

.cylVol <- function(X) {
	sum(sapply(X, function(x) { pi*x$r^2*x$length + 4/3*pi*x$r^3 }))
}

#' Plot fibre system
#'
#' Draw 3d spherocylinders
#'
#' The function requires the package \code{rgl} to be installed.
#'
#' @param S				a list of cylinders
#' @param box			simulation box
#' @param draw.axes		logical: if true, draw the axes
#' @param draw.box	    logical: if true, draw the bounding box
#' @param clipping 		logical: if true clip to the bounding box
#' @param ...			further material properties passed to 3d plotting functions
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
#
##' Cylinder vertical intersection
##' 
##' Intersect a clyinder system by a vertical section plane 
##' 
##' The function performs a vertical intersection defined by the normal vector
##' \code{n=c(0,1,0)} which depends on the main orientation axis of the
##' coordinate system and has to be parallel to this.
##'
##' @param S		 list of cylinders, see \code{\link{simCylinderSystem}}
##' @param d 	 distance of intersecting plane to the origin
##' @param n 	 normal vector of intersting plane
##' @param intern \code{intern=FALSE} (default) return all section profiles otherwise
##' 				only those which have their centers inside the intersection window
##' 
##' @return list of size, shape and angle of section profiles
#cylinderIntersection <- function(S, d, n = c(0,1,0), intern=FALSE) {
#	stopifnot(is.logical(intern))
#	if(sum(n) > 1)
#		stop("Normal vector is like c(0,1,0). ")
#	if(!(class(S) %in% c("cylinder")))
#		stop("Class must be `cylinder`.")
#	.Call(C_IntersectCylinderSystem, attr(S,"eptr"), n, d, intern, 100)	
#}
