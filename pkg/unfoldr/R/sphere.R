###############################################################################
# Author:  M. Baaske
# Date:	   2018/06/15	
# File:    sphere.R: 
# 
# Comment: simulation and intersection of sphere systems, implements EM algorithm
# 
###############################################################################


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


#' Simulation of sphere system
#'
#' The function simulates a Poisson sphere system.
#'
#' Any distribution for the radii can be specified as a name in `\code{rdist}`as long as the formal named function parameters
#' match the actual named parameters exactly as defined in the parameter list `\code{theta}`.
#'
#' The simulation box is given by a list of length three. Each element is a vector corresponding to the lower and upper points in x-, y- and z-direction in
#' the sense of a canonically oriented simulation box. If \code{box} has only one element, i.e. \code{list(c(0,1)}, the same extent
#' is used for the remaining dimensions. The argument \code{pl} denotes the print level of information during simulation and is used to set the type of the
#' return value. If \code{pl=10} only the radii are returned otherwise a full list of spheres. The argument `\code{rdist}` declares
#' the name of the (user defined) radii distribution. Setting \code{rdist="rlnorm"} leads to lognormally distributed
#' radii of the spheres. In this case an exact simulation [1] is available if \code{perfect=TRUE}. Currently available other
#' distributions are defined by `\code{rbeta}`, `\code{rgamma}`, `\code{rlnorm}` and `\code{runif}`. Use `\code{const}` for a constant
#' radius of all spheres. 
#'
#' @param theta    simulation parameters
#' @param lam      mean number of spheres per unit volume
#' @param rdist    name of radii distribution
#' @param box 	   simualtion box
#' @param dz	   distance of the intersecting x-, y-plane to the origin
#' @param n		   normal vector of intersting plane
#' @param profiles logical, \code{profiles=FALSE} (default), whether the simulated system of spheres is intersected afterwards in which case
#' 				   only sections profiles are returned 
#' @param intern   logical, \code{intern=FALSE} (default), whether to return only section profiles with centers inside the simulation window  
#' @param perfect  logical, \code{perfect=TRUE} (default), exact simulation
#' @param pl 	   integer, print level and return value definition
#' @param label    character, label passed to each simulated sphere, ´\code{N}´ (default)
#'
#' @return The function either returns a list of spheres with elements \code{id}, \code{center} and radius \code{r} of class \code{spheres}
#'  	   or a vector of radii.
#'
#' @references
#'	\itemize{
#'    \item{} {C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms.
#'             Springer, Berlin, 2002. Zbl 0990.86007}
#' 	 }
#' @examples
#'  theta <- list("meanlog"=-2.5,"sdlog"=0.2)
#'  S <- simSphereSystem(theta,lam=1000,rdist="rlnorm",pl=101)
#' 
#' @author M. Baaske
#' @rdname simSphereSystem
#' @export
simSphereSystem <- function(theta, lam, rdist, box=list(c(0,1)), dz=0, n=c(0,1,0),
							 profiles = FALSE, intern=FALSE, perfect=TRUE, pl=0, label="N")
{
	if(!is.numeric(lam) || !(lam>0) )
		stop("Expected 'lam' as non-negative numeric argument")	
	if(!is.list(theta))
		stop("Expected 'radii' as list of named arguments.")
	if(length(theta) == 0L)
		stop("Arguments for the `radii` distribution must be given.")
	if(length(box)==0 || !is.list(box))
		stop("Expected 'box' as list.")
	if(length(box)==1)
		box <- rep(box[1],3)
	if(is.null(names(box)) || !(names(box) %in% c("xrange","yrange","zrange")))
		names(box) <- c("xrange","yrange","zrange")
	if(sum(n)>1 )
	 stop("Normal vector expected as, e.g. c(0,1,0).")
	
 	cond <- list("rdist"=rdist,
				 "lam"=lam,
				 "box"=box,
				 "pl"=pl,
				 "rho"=.GlobalEnv,
			     "perfect"=as.integer(perfect),
				 "dz"=dz, "nsect"=n,
				 "intern"=as.integer(intern),
				 "label"=as.character(label))

	if(cond$rdist=="const") {
		if(length(theta)==0L)
			stop("Arguments for shape distribution must be given.")
	} else if(exists(cond$rdist, mode="function")) {
		fargs <- names(formals(cond$rdist))
		if(cond$rdist %in% c("rlnorm","rbeta","rgamma","runif"))
		  fargs <- fargs[-1]

		it <- match(names(theta),fargs)
		if(length(it)==0 || anyNA(it))
			stop(paste0("Arguments of 'theta' must match formal arguments of function: ",cond$rdist))
		
	} else
	   stop(paste("Undefined distirbution function `", cond$rdist, "` for radii."))
   	
    if(profiles) {
		.Call(C_SimulateSpheresAndIntersect, theta, cond)
	} else {
		structure(.Call(C_SphereSystem, theta, cond),
			"lam"=lam, "box" = box, "perfect"=perfect
		)	
			
	}
}

#' Sphere planar section
#'
#' Intersect sphere system
#'
#' Given a sphere system obtained from \code{\link{simSphereSystem}}
#' the function returns the section radii from intersecting all spheres
#' stored in \code{S}.
#'
#' @param S 		list of spheres of class \code{sphere}, see \code{\link{simSphereSystem}}
#' @param d 		distance of the intersecting xy-plane to the origin
#' @param intern 	logical, \code{FALSE} (default), return all planar sections otherwise
#' 					only those which have their centers inside the intersecting window
#' @param pl		print level, default \code{pl=0}
#'
#' @return 			numeric vector of circle diameters
#' 
#' @author M. Baaske
#' @rdname planarSection
#' @export
planarSection <- function(S, d, intern=FALSE, pl=0) {
   stopifnot(is.logical(intern))
   if(!(c("sphere") %in% class(S) ))
	 stop("Expected spheres as list argument.")
  sp <- .Call(C_IntersectSphereSystem,as.character(substitute(S)),c(0,0,1),d,intern,.GlobalEnv,pl)
  if(is.list(sp))								# full list of section profiles
   return (sapply(sp,function(x) 2.0*x$r))		
  else return (2*sp)							# only radii are returned
}

#' Sphere intersection
#' 
#' Intersect a sphere system
#' 
#' Given a sphere system obtained from \code{\link{simSphereSystem}}
#' the function returns the section radii from intersecting all spheres
#' stored in \code{S}.
#'
#' @param S 		list of spheres, see \code{\link{simSphereSystem}}
#' @param d 		distance of the intersecting x-, y-plane to the origin
#' @param n			normal vector of intersting plane
#' @param intern 	logical, \code{FALSE} (default), return all planar sections otherwise
#' 					only those which have their centers inside the intersecting window
#' @param pl		print level, default \code{pl=0}
#'
#' @return 			either a list of section profiles or, if \code{pl=10}, a numeric vector of radii only
#' 
#' @author M. Baaske
#' @rdname sphereIntersection
#' @export 
sphereIntersection <- function(S, d, n=c(0,1,0), intern=FALSE, pl=0) {	
	stopifnot(is.logical(intern))
	if(sum(n)>1 )
	  stop("Normal vector is like c(0,1,0). ")
    if(!(c("sphere") %in% class(S) ))
		stop(paste0("Class of argument `S` is ",class(S), " but must be `sphere`."))
	.Call(C_IntersectSphereSystem,as.character(substitute(S)), n, d, intern, .GlobalEnv, pl)	
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