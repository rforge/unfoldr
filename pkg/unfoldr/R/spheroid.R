## Comment: simulation, intersections and visualization of
## spheroid, sphere and spherocylinder 3D systems

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


#' Check intersection
#'
#' Check itersection of either spheres, spheroids or spherocylinders
#'
#' For a given list of spheres, spheroids or spherocylinders the function tests whether an object intersects
#' one of the lateral planes (bottom, top plane included) of the simulation box. The vector returned is of length
#' equal to the number of objects in \code{S} and has entries either \code{1} for an object which is intersected by
#' and \code{0} otherwise.
#'
#' @param S 	list of spheres, spheroids or spherocylinders, see \code{\link{simPoissonSystem}}
#'  
#' @return 		binary integer vector of length equal to the length of \code{S}
#'  
#' @author M. Baaske
#' 
#' @rdname updateIntersections
#' @export
updateIntersections <- function(S) {
	.Call(C_UpdateIntersections, as.character(substitute(S)), .GlobalEnv)
}

#' Construct section profiles
#'
#' Set up section profiles of spheroids for unfolding
#'
#' The function aggregates the necessary information for trivariate unfolding of spheroids' joint size-shape orientation distribution
#' of type either "\code{prolate}" or "\code{oblate}". The argument \code{size} is a numeric matrix of semi-axis lengths where the first column
#' corresponds to the major semi-axis and the second one to minor semi-axis. The orientation of an ellipse is assumed to be measured as the angle
#' between its major axis and vertical axis of the coordinate system in the intersection plane ('z' axis in 3D). For values in \eqn{[0,2\pi]} these
#' angles are automatically transformed to \eqn{[0,\pi/2]} as required by the unfolding procedure.
#'
#' @param size	  matrix of lengths of the semi-axes
#' @param alpha   angle of section profiles in the plane (see details)
#' @param type    name of the spheroid type, either "\code{prolate}" or "\code{oblate}" from which the
#'                section profiles are assumed to come from
#'
#' @return		The function returns a list which consists of either the longer or shorter
#' 				semi-axis length named \code{A} of section profiles corresponding to the type of spheroids used before and whose joint
#' 				joint distribution is to be estimated (by unfolding), the shape factor \code{S} of both semi-axes as the shape factor
#' 				between \eqn{(0,1]} and the orientation	angle \code{alpha}, either of class "\code{prolate}" or "\code{oblate}".
#'
#' @examples
#'  # load data set
#'  data(data15p)
#'  
#'  # matrix of semi-axes lengths (major,minor)
#'  AC <- data.matrix(data15p[c("A","C")])/1000	
#' 
#'  # selecting the minor semi-axis for prolate type of spheroids:
#'  # independent of nomenclature (always named \code{A})
#'  sp <- sectionProfiles(AC,unlist(data15p["alpha"]))
#' 
#'  summary(sp$A)			# here minor semi-axis because of prolate
#'  summary(sp$S)			# shape factor
#'  summary(sp$alpha)		# angle assumed to be w.r.t. (vertical) 'z' axis 
#'  
#' @author M. Baaske
#' @rdname sectionProfiles
#' @export
sectionProfiles <- function(size,alpha,type=c("prolate","oblate")) {
	type <- match.arg(type)
	stopifnot(is.matrix(size))
	if(anyNA(size) || any(size<0))
		stop("'size' must have non-negative values.")
	if(anyNA(alpha) || !is.numeric(alpha) || any(alpha<0))
		stop(paste("'alpha' must have non-negative values."))
	if(max(alpha) > pi/2)
	 alpha <- try(sapply(alpha,.getAngle),silent=TRUE)
    if(inherits(alpha,"try-error"))
	 warning("Could not compute 'alpha'. See its returned list element.")	
    
    structure(list("A"=if(type=="prolate") size[,2] else size[,1],
				   "S"=size[,2]/size[,1],
				   "alpha"=alpha),							
		   class=type)
}

#' Poisson germ-grain process
#'
#' Simulation of Poisson germ-grain processes with either spheres, spheroids or spherocylinders as grains
#'
#' The function can simulate a Poisson germ-grain process according to the parameter \code{theta} within a predefined (3D) box.
#' The positions of the germs follow a uniform distribution according to a Poisson process with mean intensity parameter \code{lam}.
#' The function either randomly generates \code{type="prolate"} or \code{type="oblate"} spheroids, spheres or spherocylinders.
#' The argument \code{size} sets the name of the distribution function for the size/length of the objects, i.e. the major semi-axis
#' lengths in case of spheroids, radii for spheres or the lengths of the main axis of rotation for spherocylinders including the end caps. 
#' 
#' The following direction (orientation) distributions of the spheroids' major-axis, respectively, cylinders' main axis are available:
#' a uniform distribution ("\code{runifdir}"), distribution ("\code{rbetaiso}") and the "\emph{von Mises-Fisher}" ("\code{rvMisesFisher}")
#' distribution. The two last ones depend on the concentration parameter \code{kappa} which is set as part of the parameter list \code{theta}, see examples below.
#' The direction distributions generate random spherical coordinates \eqn{(\vartheta, \phi)} w.r.t. a fixed main orientation axis \code{mu}
#' with polar angle \eqn{\vartheta\in[0,\pi/2)} and azimuthal angle \eqn{\phi\in[0,2\pi)}. The simulations are always performed within a bounding 3D box which consists of a list specifying the ranges of each dimension corresponding to the lower and upper limits of the box in each direction. If the
#' argument \code{box} contains only a single range, i.e. \code{box=list(c(0,1))}, this limit isassumed for the remaining dimensions which is then simply replicated.
#' The optional argument \code{rjoint} defines a (joint) distribution function which can be any function provided by the user in order to generate
#' the required distributional parameters for the spheroids or spherocylinders. For an in-depth example of usage please see the workflow in 'simSpheroids.R'
#' and 'simCylinders.R'. 
#' 
#' In addition, the function supports an exact simulation type [2] of the grains. In case of spheroids and spherocylinders setting \code{size="rbinorm"}
#' declares a bivariate size-shape distribution for which the exact simulation is available. More specifically, for a bivariate normal random vector \eqn{[X,Y]}
#' with correlation parameter \eqn{\rho}, the length of the major semi-axis of a spheroid is given by \eqn{a=exp(x)} with a (logit-transformed) shape parameter
#' as \eqn{s=1/(1+exp(-y))} and thus a scaled minor semi-axis length \eqn{c=a*s}. The modification leads to a log-normally distributed length of the
#' major semi-axis. Consequently, in case of spherocylinders, the log-normally distributed length is \eqn{len=h+2*r} where \code{h} is the height and
#' \eqn{r=len/2*s} the radius. The main direction \code{u} of a spheroid or spherocylinder is determined by the major axis independent of size and shape.
#' Further, the following univariate distributions of the major semi-axis \code{a}, respectively, length \code{len} and shape \code{s} are available:
#' '\code{rbeta}', '\code{rgamma}', '\code{rlnorm}' and '\code{runif}'. One can also use '\code{const}' for simulations with constant lengths or shapes.
#' Note that only simulations with size distributions '\code{rbinorm}' or '\code{rlnorm}' can use the exact type of simulation.
#'
#' For spheres any distribution of the radii can be specified as a name of a user-defined function in the argument \code{size} as long as the formal named
#' function parameters match the actual names of the parameters exactly as defined in \code{theta}. Besides this, all other distributions given above are
#' also available. Using '\code{const}' simulates spheres of constant radii. 
#'  
#' The argument \code{pl>=0} denotes both the print level of intermediate output and by the same time the type of the return value. If \code{pl=10},
#' then an abbreviated list of spheroids or spheres is returned to speed up computation. Note that, the current implementation does not include routines
#' for unfolding the joint size-shape-orientation distribution of spherocylinders so far. 
#' 
#'
#' @param theta 	list of simulation parameters which must consist of elements: \code{size}, \code{shape} and \code{orientation}
#' @param lam   	mean number of objects per unit volume
#' @param size  	name of the size distribution function 
#' @param shape 	name of the shape distribution function  
#' @param orientation name of direction distribution function
#' @param type     type of grain, either "\code{prolate}" or "\code{oblate}", "\code{spheres}", "\code{cylinders}"
#' @param rjoint   user-defined function, which specifies the (joint) distribution of the size, shape and orientation 
#' @param box	   simulation box
#' @param mu	   main orientation axis, \code{mu=c(0,0,1)} (default)
#' @param dz	   distance of the intersecting plane to the origin
#' @param n		   normal vector defining the intersecting plane
#' @param intersect options for type of return values: "\code{full}" for the simulated system together with section profiles as two lists named \code{S} and \code{sp} respectively,
#'                  choose "\code{only}" for section profiles only, "\code{original}" for the 3D system only and "\code{digi}" for a (binary) integer matrix \code{W} as a discretized
#' 				    version of section profiles whose resolution depends on the chosen lattice constant \code{delta}
#' @param delta	   lattice constant for discretization, set to \code{0.01} (default)
#' @param intern   logical, \code{FALSE} (default), whether to return only section profiles with centers inside the simulation window
#' @param perfect  logical, \code{FALSE} (default), whether to simulate exactly (also called perfect)
#' @param pl  	   integer, print level and return value type, see details
#' @param label    character, a label set to each generated object, set to '\code{N}' (default)
#'  
#' @return 		   list of 3D objects depending on the chosen return type defined by the argument \code{intersect}
#'
#' @examples
#'  # intensity parameter
#'  lam <- 100
#' 
#'  # simulation bounding box
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#' 
#'  # log normal size distribution with a constant shape factor and
#'  # concentration parameter (\code{kappa=1}) for the orientation, see reference [1] 
#'  theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5),
#'                "shape"=list("s"=0.5),
#'                "orientation"=list("kappa"=1))
#' 
#'  S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="oblate",pl=1) 
#'  length(S)
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
								type=c("prolate","oblate","spheres","cylinders"), rjoint=NULL,
								  box=list(c(0,1)), mu=c(0,0,1), dz=0, n=c(0,1,0),
								  intersect=c("original","full","only","digi"), 
								    delta=0.01, intern=FALSE, perfect=FALSE, pl=0, label="N")
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
	if(dz < min( box[[which(n == 1)]]) || dz > max(box[[which(n == 1)]]))
	 stop("The plane has to intersect the box.")
	
	# digitization (resolution) factor
	stopifnot(is.numeric(delta))
	
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
			 "dz"=dz, "nsect"=n, "delta"=delta,"intern"=as.integer(intern),
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
			# set defaults
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
				stop(paste0("Arguments of 'shape' must match formal arguments of function ",shape,sep=""))
		} else {
			stop(paste0("Undefined `", shape, "` distribution function."))
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
				 stop(paste0("Arguments of 'size' must match formal arguments of function ",size,sep=""))
			 		
		 } else if(size == "const"){
			 # set defaults
			 if(!is.list(theta$size) || length(theta$size) == 0L)
			   theta$size <- list(1,1)
		 } else {
			stop(paste0("Undefined `", size, "` distribution function."))
		 }
		 
		 it <- pmatch(orientation,c("runifdir","rbetaiso","rvMisesFisher","const"))
		 if(is.na(it) && !exists(orientation, mode="function"))
			 stop("Undefined distribution function for the orientation/direction.")	 
		 
		 # set defaults
		 if(!is.list(theta$orientation) || length(theta$orientation) == 0L)
			 theta$orientation <- list("kappa"=1)
		 
		 list("type"=type, "lam"=lam,
			  "rdist"=list("size"=size, "shape"=shape, "orientation"=orientation),
			  "box"=box,"pl"=pl,"mu"=mu,"rho"=.GlobalEnv, "dz"=dz, "nsect"=n, "delta"=delta,
			  "intern"=as.integer(intern),"label"=label, "perfect"=as.integer(perfect),
			  "intersect"=intersect)
	}
		 
	.Call(C_PoissonSystem, theta, cond)	
}

#' Coefficients for trivariate unfolding
#'
#' Compute coefficients of the discretized integral equation for unfolding
#'
#' In order to apply the Expectation Maximization (EM) algorithm to the stereological unfolding procedure for the joint
#' size-shape-orientation distribution in 3D one first has to compute the coefficients of a discretized integral equation.
#' This step is the most time consuming part of unfolding and therefore has been separated in its own function and can be called
#' separately if needed. The number of bin classes for the size, shape and orientation do not need to be the same, but the
#' given class limits are also used for binning the spatial values. One can define the number of cpu cores by the global option
#' \code{par.unfoldr} or passing the number of cores \code{nCores} directly to the function.
#'
#' @param   breaks  list of bin vectors, see \code{\link{setbreaks}}
#' @param   stype   type of spheroid, either "\code{prolate}" or "\code{oblate}"
#' @param   check   logical, whether to run some input checks
#' @param   nCores  number of cpu cores used to compute the coefficients
#' 
#' @return  numeric 6D array of coefficients
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
			stop(paste0("Breaks vector 'angle' must have values between zero and ",pi/2))
		if(min(breaks$shape)<0 || max(breaks$shape)>1)
			stop("Breaks vector 'shape' must have values between 0 and 1.")
	}

	.Call(C_CoefficientMatrixSpheroids,
			breaks$size,breaks$angle,breaks$shape,
			breaks$size,breaks$angle,breaks$shape, list(stype,nCores))

}

#' Vertical sections
#'
#' Compute vertical section profiles of a spheroid system
#'
#' The function intersects a spheroid system by a plane defined by the normal vector \code{n} either
#' equal to \code{c(0,1,0)} (default) or \code{c(1,0,0)}, which is called a vertical section. Depending on
#' the type of spheroid (either "\code{prolate} or "\code{oblate}") the returned semi-axis lengths are those
#' corresponding to the minor semi-axis or, respectively, major semi-axis in the way these are required for unfolding. 
#'
#' @param S		 list of spheroids, see \code{\link{simPoissonSystem}}
#' @param d 	 distance of the intersecting plane from the origin of the box
#' @param n 	 normal vector which defines the interecting vertical plane
#' @param intern logical, \code{FALSE} (default), return all section profiles otherwise
#' 				 only those which have their centers inside the correspondig intersection window
#' 
#' @return 	 	 list of sizes \code{A}, shape factors \code{S} and (vertical) angles \code{alpha}
#'               of section profiles in the plane w.r.t the 'z' axis between \eqn{[0,\pi/2]}.
#' 
#' 
#' @examples  
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#'  
#'  # (exact) bivariate size-shape (isotropic) orientation distribution (spheroids)
#'  theta <- list("size"=list("mx"=-2.5,"my"=0.5, "sdx"=0.35,"sdy"=0.25,"rho"=0.15),
#' 		"orientation"=list("kappa"=1))
#' 
#'  S <- simPoissonSystem(theta,lam=100,size="rbinorm",box=box,
#'   type="prolate",perfect=TRUE,pl=1)
#' 
#'  sp <- verticalSection(S,d=2.5,n=c(0,1,0),intern=TRUE)
#'  summary(sp$alpha)
#'  
#' @author M. Baaske
#' @rdname verticalSection 
#' @export
verticalSection <- function(S,d,n=c(0,1,0),intern=FALSE) {
	if(!(class(S) %in% c("prolate","oblate")))
	 stop("Only spheroids of class 'prolate' or 'oblate' can be used.")
 	stopifnot(is.numeric(d))	
 	stopifnot(is.logical(intern))
	if(sum(n) != 1 || which(n == 1) > 2L)
	  stop("Normal vector is one of c(0,1,0) or c(1,0,0).")
	
    ss <- .Call(C_IntersectPoissonSystem,
		 		 as.character(substitute(S)),
		  		 list("nsect"=n,"dz"=d,"intern"=intern,"pl"=10),
		  		 .GlobalEnv)
  	
	A <- if(class(S)=="prolate")
			sapply(ss,"[[",2)
		 else sapply(ss,"[[",1)
	
    ## convert angle 'alpha' in the intersecting plane
	## which is always between [0,pi/2] and w.r.t 'z' axis
	alpha <- sapply(ss,"[[",4)
    alpha <- try(0.5*pi-sapply(alpha,.getAngle),silent=TRUE)		# alpha in [0,pi/2]
	if(inherits(alpha,"try-error") || !is.numeric(alpha) || anyNA(alpha) )
	   warning("Could not compute angle 'alpha'. See the returned list element.")
    	
	structure(list("A"=A,
				   "S"=sapply(ss,"[[",3),
				   "alpha"=alpha),	# w.r.t z axis in 3D
		   "win"=attr(ss,"win"),"plane"=attr(ss,"plane"),
		 class=class(S) )
}

#' Intersection in 3D
#'
#' Intersect a system of spheres, spheroids or spherocylinders by a plane
#'
#' The function intersects a given (Poisson) system made of spheres, spheroids or cylinders as grains by a plane defined by a
#' normal vector, e.g. \code{n=c(0,1,0)}, perpendicular to one of the bounding planes of the simulation box. 
#' For a print level \code{pl>=0} some verbose output is given. Also it sets the type of return value. In case of spheroid
#' intersections, setting \code{pl=10}, leads to a short version of the full specification return list of section profiles
#' with elements named \code{A} (major semi-axis), \code{C} (minor semi-axis), \code{S} (the shape factor as the ratio of these two)
#' and \code{phi} as the angle in the intersecting plane between \eqn{[0,2\pi]} w.r.t. the 'x' axis. Otherwise additional components
#' are returned such as the ellipse`s rotation matrix, also named \code{A}, the center point \code{center} and a constant \code{type=10}
#' (defining the object of full ellipses among other types of intersection objects) of the section profiles. For sphere intersections only
#' a numeric vector of disc radii are returned as a short version of return values.     
#'
#' @param S		 list of spheres, spheroids or spherocylinders, see \code{\link{simPoissonSystem}}
#' @param d 	 distance of the the box-aligned intersecting plane from the origin
#' @param n 	 normal vector which defines the intersecting plane
#' @param intern logical, \code{FALSE} (default), return all section profiles otherwise
#' 				 only those which have their centers inside the corresponding intersection window
#'               (if the intersected system had been simulated using exact simulation this makes sense)
#' @param pl	 integer, \code{pl=0} (default), for no verbose output and otherwise a full specification list of
#'				 section profiles in case of sphere and spheroid intersections
#' 				
#' @return 		For spheroid intersections the function returns a list of size, shape and angle of section profiles or a short version
#' 				of it; for sphere intersections either radii or lists containing the centers of discs and the object number.
#' 
#' 
#' @examples
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#'  
#'  # constant size-shape orientation distribution (spheroids)
#'  theta <- list("size"=list(0.1),"shape"=list(0.5), "orientation"=list("kappa"=10))
#' 
#'  S <- simPoissonSystem(theta,lam=100,box=box,type="prolate",pl=1)
#'  
#'  # return short version of section profiles
#'  sp <- intersectSystem(S, 2.5, pl=10)		
#' 
#' @author M. Baaske
#' @rdname intersectSystem
#' @export
intersectSystem <- function(S, d, n=c(0,1,0), intern=FALSE, pl=0) {
	if(!(class(S) %in% c("prolate","oblate","spheres","cylinders")))
	  stop("Unknown object class.")
    stopifnot(is.numeric(d))
	stopifnot(is.numeric(pl))
	stopifnot(is.logical(intern))
	if(sum(n)>1)
	  stop("Normal vector is like c(0,1,0). ")
    box <- attr(S,"box")
	if(is.null(box))
	 stop("'box' must be given as an attribute.")
    if(d < min( box[[which(n == 1)]]) || d > max(box[[which(n == 1)]]))
	  stop("The plane has to intersect the box.")
	
	.Call(C_IntersectPoissonSystem,
			as.character(substitute(S)),
			list("nsect"=n,"dz"=d,"intern"=as.integer(intern),"pl"=as.integer(pl)),
		    .GlobalEnv)	
}

#' Spheroid system 3D
#'
#' Draw a spheroid system in 3D
#'
#' The function requires the package '\code{rgl}' to draw spheroids. For an example
#' please see the file 'simSpheroids.R'.
#'
#' @param S				list of spheroids, see \code{\link{simPoissonSystem}}
#' @param box			simulation box
#' @param draw.axes		logical, whether to show the axes
#' @param draw.box	    logical, whether to show the bounding box
#' @param draw.bg	    logical, whether to show a gray background box correpsonding to the simulation box
#' @param bg.col 		background color used showw the box background
#' @param clipping 		logical, whether to clip all lateral planes of the simulation box
#' @param ...			further material properties passed to 3D plotting functions
#' 
#' @return NULL
#' @author M. Baaske
#' @rdname spheroids3d
#' @export
spheroids3d <- function(S, box, draw.axes=FALSE, draw.box=TRUE, draw.bg=TRUE,
		           bg.col="white", clipping=FALSE, ...)
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


#' Cylinder system 3D
#'
#' Draw spherocylinders in 3D
#'
#' The function requires the package '\code{rgl}' to draw spherocylinders into a 3D '\code{rgl}' image.
#' For an example please see the file 'simCylinders.R'.
#'
#' @param S				list of cylinders, see \code{\link{simPoissonSystem}}
#' @param box			simulation box
#' @param draw.axes		logical, whether to show the axes
#' @param draw.box	    logical, whether to show the bounding box
#' @param clipping 		logical, whether to clip all lateral planes of the simulation box
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
	ok <- sapply(S, function(x) (!is.null(x$h) && x$h > 0.0) ) 		# only true cylinders not 'spheres'
	
	cyls <- lapply(S[ok], function(x) { cylinder(x$center, x$r, x$h, x$rotM, x$u) })
	rgl::shapelist3d(cyls,...)
	
	if("col" %in% names(args)) {
		cols <- rep(rep(args$col,length.out=length(S)),each=2)
		args$col <- NULL
	} else cols <- "black"
	
	# spheres
	is.sphere <- sapply(S,function(x) !(is.null(x$h) || x$h == 0.0) )
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


#' Spheroid intersections 3D
#'
#' Draw section profiles of spheroids in 3D
#'
#' The function requires the package '\code{rgl}' for drawing section profiles (ellipses)
#' into a 3D '\code{rgl}' image. For a full example please see the file 'simSpheroids.R'.
#'
#' @param E				a list of spheroid intersections, see \code{\link{intersectSystem}}
#' @param n 			the normal vector of the intersecting plane
#' @param np			number of points for a polygon approximation of the ellipses
#' 
#' @return NULL
#' 
#' @author M. Baaske
#' @rdname drawSpheroidIntersection
#' @export
drawSpheroidIntersection <- function(E, n=c(0,1,0), np=25) {
	ind <- which(n == 0)
	
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

#' Sphere planar sections
#'
#' Planar intersection of spheres
#'
#' The function computes the planar intersection of a sphere system, i.e. an intersection with the plane whose normal
#' vector is given by \code{c(0,0,1)} and returns the diameters of the resulting discs.
#'
#' @param S 		list of spheres of class \code{sphere}, see \code{\link{simPoissonSystem}}
#' @param d 		distance of the (planar) xy-plane to the origin
#' @param intern 	logical, \code{FALSE} (default), whether to return all discs or
#' 					only those which have their centers inside the intersecting window
#' @param pl		print level, \code{pl>0} for some verbose output
#'
#' @return 			numeric vector of disc diameters
#' 
#' @examples
#'  lam <- 100
#'  
#'  # parameter rlnorm distribution (radii)
#'  theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5))
#'  
#'  # simulation bounding box
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#'  
#'  # simulate only 3D system
#'  S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="spheres",
#'    intersect="original", pl=1)
#'  
#'  # return only objects whose centers are within
#'  # the intersection window
#'  sp <- planarSection(S,d=2.5,intern=TRUE,pl=1)
#'  
#'  # histogram of diameters
#'  hist(sp)
#'  summary(sp)
#'  
#'  # distribution of radii
#'  mean(log(sp/2))
#'  sd(log(sp/2))
#' 
#' 
#' @author M. Baaske
#' @rdname planarSection
#' @export
planarSection <- function(S, d, intern=FALSE, pl=0) {
	if(!(c("spheres") %in% class(S) ))
	 stop("Expected spheres as a list.")
	if(!is.numeric(d) || d<0)
	 stop("'d' must be numeric/positive.")
	stopifnot(is.logical(intern))
	
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
#' The function provides basic binning (grouping) of numeric values into
#' classes defined by the breaks vector \code{bin}. The values are binned
#' according to \eqn{bin[i[j]]<x[j]\leq bin[i[j]+1]} for intervals \eqn{i=1,...,N-1}
#' and \code{length(bin)=N} of values \eqn{x[j]}, \eqn{j=1,...,|x|}. If \eqn{x[j] > bin[N]} or \eqn{x[j] < bin[1]} then \eqn{x[j]}
#' is not counted at all.
#'
#' @param x   	 numeric values to be binned
#' @param bin    non-decreasingly sorted breaks vector
#' @param na.rm  logical, default \code{FALSE}, whether to remove missing values, including \code{NaN} in \code{x}
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

#' Coefficients for Expectation Maximization algorithm
#'
#' Matrix of coefficients for Wicksell's corpuscle problem
#'
#' The function computes the matrix of coefficients for solving the discretized integral equation of
#' the Wicksell's corpuscle problem.
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

#' Expectation Maximization (EM) algorithm
#'
#' Estimation of the empirical sphere diameter distribution
#'
#' The function performs the EM algorithm, see reference below.
#'
#' @param y   		vector of observed absolute frequencies of circle diameters
#' @param bin 		non-decreasing vector of class limits
#' @param maxIt		maximum number of iterations to be used
#' @return 			vector of count data of absolute frequenties of sphere diameters
#'
#' @example inst/examples/simSpheres.R
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

#' Digitization
#'
#' Digitize 2D section profiles  
#'
#' A list of section profiles, e.g. discs, ellipses or segments from intersected spherocylinders, is digitized according
#' to a given resolution according to the lattice constant \code{delta} such that the result is an integer matrix 
#' which can be interpretated as a binary image. An intersection window can be either provided by the user w.r.t. \eqn{[l,u]^2}
#' where \code{l,u}	are lower, respectively, upper bounds of corresponding to a simulation box, or, it is taken from the section
#' profiles \code{sp} which stores the intersection window as an attribute.
#'
#' @param sp   		list of section profiles, see \code{\link{intersectSystem}}
#' @param delta		lattice constant for discretization of section profiles
#' @param win		list of length two, the intersection window, default \code{NULL} 
#' 
#' @examples
#'  # simulation box
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#'  # (exact) bivariate size-shape (isotropic) orientation distribution (spheroids)
#'  theta <- list("size"=list("mx"=-2.5,"my"=0.5, "sdx"=0.35,"sdy"=0.25,"rho"=0.15),
#'                "orientation"=list("kappa"=1))
#' 
#'  # return only 3D system
#'  S <- simPoissonSystem(theta,lam=100,size="rbinorm",box=box,type="prolate",
#'        intersect="original",n=c(0,1,0),mu=c(0,0,1),perfect=TRUE,pl=1)
#' 
#'  # vertical intersection w.r.t. 'mu' (z axis, see above)
#'  sp <- intersectSystem(S, 2.5)
#' 
#'  # show intersecting window
#'  win <- attr(sp,"win")
#' 
#'  # digitize (could also pass some 'win' as an argument) 	
#'  W <- digitizeProfiles(sp, delta=0.01, win = NULL)  
#'  image(1:nrow(W),1:ncol(W),W,col=gray(1:0))
#' 
#' @return 			binary (integer) matrix as an image 
#'
#' @author M. Baaske
#' @rdname digitizeProfiles
#' @export 
digitizeProfiles <- function(sp, delta = 0.01, win = NULL ) {
	stopifnot(is.numeric(delta))
	if(!is.null(win)){
		if(!is.list(win) || length(win) == 0L )
			stop("Expected argument 'win' as list type.")
		if(length(win) == 1L)
			win <- rep(win[1],2)
		else if(length(win) != 2L)
			stop("Digitization window 'win' must be two-dimensional.")		
	}
	if(delta < 0 || delta > 1.0)
	  stop("'delta' should be positive and equal or less than 1.")
  
	.Call(C_DigitizeProfiles,as.character(substitute(sp)), delta, win, .GlobalEnv)
}