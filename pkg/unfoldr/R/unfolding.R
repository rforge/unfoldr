## Comment: 
## Functions to estimate the joint size-sjape-orientation distribution
## of prolate or oblate spheroids (ellipsoids of revolution),
## visualization of the trivariate 'unfolded' histogram of size, shape
## and orientation, implements the EM algorithm for binned data for 
## spheres and spheroids (not yet for spherocylinders)


#' Trivariate stereological unfolding
#'
#' Estimate the joint size-shape-orientation distribution of spheroids
#'
#' Given an array of coefficients \code{P}, see \code{\link{coefficientMatrixSpheroids}} and an input histogram
#' \code{F} of measured planar characteristics of section profiles, the function estimates the spatial joint
#' size-shape-orientation distribution of the corresponding spheroids in 3D by a discretized version of the
#' \emph{Expectation Maximization} (EM) algorithm. A number of cpu cores can be set by the option '\code{par.unfoldr}'
#' for parallel computations. The function is also internally called by \code{\link{unfold}} in case of spheroids.
#'
#' @param P 		coefficient array
#' @param F 		input histogram
#' @param maxIt 	maximum number of EM iterations
#' @param nCores 	number of cpu cores to be used
#' 
#' @return trivariate histogram
#'
#' @example inst/examples/unfold.R
#'
#' @references
#' 	Bene\eqn{\check{\textrm{s}}}, V. and Rataj, J. Stochastic Geometry: Selected Topics Kluwer Academic Publishers, Boston, 2004
#' 
#' @author M. Baaske
#' @rdname em.spheroids
#' @export
em.spheroids <- function(P,F,maxIt,nCores=getOption("par.unfoldr",2L)) {
	.Call(C_EMS,P,F,list("maxSteps"=maxIt,"nCores"=nCores))
}

#' Stereological unfolding
#'
#' Unfolding the (joint) distribution of planar parameters
#'
#' This is a S3 method for either trivariate stereological unfolding or estimation of the 3D diameter distribution
#' of spheres which is better known as the \emph{Wicksell's corpuscle problem}. The function aggregates all intermediate
#' computations required for the unfolding procedure given the data in the prescribed format, see reference of functions below,
#' and returning the characteristics as count data in form of a \emph{trivariate} histogram. The section profile objects \code{sp},
#' see \code{\link{sectionProfiles}}, are either of class \code{prolate} or \code{oblate} for the reconstruction of the corresponding
#' spheroids or, respectively, spheres. The result of the latter is simply a numeric vector of circle diameters. The number of bin
#' classes for discretization of the underlying integral equations which must be solved is set by the argument \code{nclass}.
#' In case of Wicksell's corpuscle problem (spheres as grains) this is simply a scalar value denoting the number of bins for the diameter.
#' For spheroids it refers to a vector of length three defined in the order of the number of size, angle and shape class limits which are used.
#' If \code{sp} is a numeric vector (such as for the estimation of the 3D diameter distribution from a 2D section of spheres) the function calls
#' the EM algorithm as described in [3].
#' The return value of the function is an object of class "\code{unfold}" with elements as follows
#' \itemize{
#' 	\item{N_A}{ (trivariate) histogram of section profile parameters}
#'  \item{N_V}{ (trivariate) histogram of reconstructed parameters}
#'  \item{P}{ array of coefficients}
#'  \item{breaks}{ list of class limits for binning the parameter values}
#' }
#'
#' @param sp 	  section profiles, see \code{\link{sectionProfiles}}
#' @param nclass  number of classes, see details 
#' @param maxIt   maximum number of EM iterations
#' @param nCores  number of cpu cores
#' @param ...	  optional arguments passed to \code{\link{setbreaks}}
#' 
#' @return        object of class "\code{unfold}", see details
#'
#' @seealso \code{\link{setbreaks}}, \code{\link{binning3d}} 
#' 
#' @examples
#'  lam <- 100
#'  # parameter rlnorm distribution (radii)
#'  theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5))
#' 
#'  # simulation bounding box
#'  box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))
#'  # simulate only 3D system
#'  S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="spheres",
#'    perfect=TRUE, pl=1)
#' 
#'  # intersect
#'  sp <- planarSection(S,d=2.5,intern=TRUE,pl=1)
#' 
#'  # unfolding
#'  ret <- unfold(sp,nclass=25)
#'  cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")
#' 
#' @author M. Baaske
#' @rdname unfold
#' @export
unfold <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",2L),...) UseMethod("unfold", sp)

#' @method unfold oblate 
#' @export
unfold.oblate <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",2L),...) {
	unfold.prolate(sp,nclass,maxIt,nCores,...)
}

#' @method unfold prolate 
#' @export
unfold.prolate <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",2L),...) {
	if(length(nclass)!=3)
	  stop("Lenght of 'dims' not equals 3.")
  	if(any(!(c("A","S","alpha") %in%  names(sp))))
		stop("Missing named arguments in 'X'")

	breaks <- setbreaks(nclass=nclass,maxSize=max(sp$A),...)
	N_A <- binning3d(sp$A,sp$alpha,sp$S,breaks)
	P <- coefficientMatrixSpheroids(breaks,class(sp),TRUE,nCores)

	N_V <- em.spheroids(P,N_A,maxIt,nCores)
	structure(list("N_A"=N_A,"N_V"=N_V,"P"=P,"breaks"=breaks),
			class=c("unfold",class(sp)))
}

#' @method unfold numeric 
#' @export
unfold.numeric <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",2L),...) {
	## The function calls the saltykov algorithm for spheres
	## (Wicksell's corpuscle problem)
	
	if(anyNA(sp))
		stop("Vector of radii 'sp' has NAs.")
	if(!is.numeric(nclass) || length(nclass) != 1L)
		stop("Expected numeric value 'nclass' as number of classes.")
	breaks <- seq(0,max(sp), length.out=nclass)

	# Input histogram	
	y <- binning1d(sp,breaks)

	n <- length(y)
	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(breaks),
			p=as.vector(matrix(0,n,n)))$p

	theta <- y+1.0e-6
	theta <- .C(C_em_saltykov,as.integer(n), as.integer(maxIt),
				as.numeric(p),as.numeric(y),theta=as.numeric(theta))$theta

	structure(list("N_A"=y,"N_V"=theta,"P"=p,"breaks"=breaks),
			class=c("unfold",class(sp)))
}

#' Histogram data
#'
#' Count data of size, shape and orientation
#'
#' For each value of planar or spatial measured quantities \code{size}, \code{shape} and \code{orientation}
#' the function counts the number of observations falling into each class (bin). The list \code{breaks} can be
#' obtained by the function \code{\link{setbreaks}}. If \code{check=TRUE}, some checks on \code{breaks} are done.
#' For an example please see the file 'almmc.R'.
#'
#' @param size  vector of sizes
#' @param angle vector of angles
#' @param shape	vector of shape factors
#' @param breaks list of bin vectors
#' @param check logical, default is \code{TRUE}
#' @param na.rm logical, whether \code{NAs} should be removed, default \code{TRUE}
#' 
#' @return 3D array of binned count data
#' 
#' @author M. Baaske
#' @rdname binning3d
#' @export
binning3d <- function(size,angle,shape,breaks,check=TRUE,na.rm = TRUE) {
	if(na.rm) {
		size <- na.omit(size)
		angle <- na.omit(angle)
		shape <- na.omit(shape)
	} else if(anyNA(c(size,angle,shape)))
		stop("Data vectors contain missing values or NAs.")

	if(!is.list(breaks) || length(breaks)!=3)
		stop("Expected 'breaks' argument as list of break vectors of length 3 for 'size', 'angle' and 'shape' data.")

	it <- match(names(breaks), c("size","angle","shape"))
	if (anyNA(it))
		stop("Expected 'breaks' as named list of: 'size','angle','shape' ")

	if (is.unsorted(breaks$size) || is.unsorted(breaks$angle) || is.unsorted(breaks$shape))
		stop("'breaks' list must contain non-decreasingly sorted classes")

	if(check) {
		if(any(breaks$size<0))
			stop("Breaks vector 'size' must have non-negative values.")
		if(min(breaks$angle)<0 || max(breaks$angle)>pi/2)
			stop(paste("Breaks vector 'angle' must have values between zero and ",quote(pi/2),sep=""))
		if(min(breaks$shape)<0 || max(breaks$shape)>1)
			stop("Breaks vector 'shape' must have values between 0 and 1.")
	}

	# return object of class 'triHist'
	.Call(C_Binning3d,size,angle,shape,breaks$size,breaks$angle,breaks$shape)
}


#' Break vectors
#'
#' Construct class limits vectors
#'
#' The function constructs the class limits for the size, shape and orientation parameters.
#' One can either define "\code{linear}" class limits of the sizes as \eqn{a_i=i\delta, \delta=maxSize/M } or
#' using exponentially increasing limits like \eqn{base^i, i=1,\dots,M}.
#' The orientation classes are defined by \eqn{\theta_j=j\omega, \omega=\pi/(2N), j=1,\dots,N} in the range
#' \eqn{[0,\pi/2]}, where \eqn{M,N} are the number of size classes and, respectively, the number of orientation classes.
#' The argument \code{base} must not be \code{NULL} if \code{sizeType} equals "\code{exp}".
#'
#' @param nclass 	number of classes
#' @param maxSize 	maximum of \code{size} values
#' @param base 		constant for size class construction
#' @param kap  	    constant for shape class construction
#' @param sizeType 	either \code{linear} or \code{exp}, default is \code{linear}
#' 
#' @examples 
#'   setbreaks(c(8,5,7),0.935,base=0.5,kap=1.25,sizeType="exp")
#'   
#' @return 			list of class limits vectors
#' 
#' @author M. Baaske
#' @rdname setbreaks
#' @export
setbreaks <- function(nclass,maxSize,base=NULL,kap=1,sizeType=c("linear","exp")) {
	if(length(nclass)==0 || any(nclass==0))
		stop("Number of bin classes 'nclass' must be greater than zero.")
	type <- match.arg(sizeType)

	binA <- switch(type,
		linear=seq(from=0,to=max(maxSize),by=max(maxSize)/nclass[1]),
		exp={
			if(is.null(base))
			 stop("Argument 'base' must be non NULL if 'sizeType' equals 'exp'.")
			ll <- unlist(sapply(1:nclass[1],function(i) base^i))
			if(base<1) ll <- sort(ll)
		  	ll
		})
	structure(list("size"=binA,
		"angle"=seq(from=0,to=pi/2,by=pi/(2*nclass[2])),
		"shape"=unlist(lapply(0:nclass[3],function(i) (i/nclass[3])^kap))))
}

#' 3D characteristics of spheroids
#'
#' Extract size, shape and orientation angle
#'
#' The function extracts the characteristics of a 3D spheroid system, either oblate or prolate spheroids,
#' and returns a list of the following elements: major semi-axis length '\code{a}', shape factor '\code{s}'
#' and polar angle '\code{Theta}'. For a full example please see the file 'unfold.R'.
#'
#' @param   S list of spheroids
#' 
#' @return  a list of spheroids` 3D characteristics
#'  
#' @author M. Baaske
#' @rdname parameters3d
#' @export
parameters3d <- function(S) {
	stopifnot(class(S) %in% c("oblate","prolate"))
	idx <- if(class(S)=="prolate") c(1,3) else c(3,1)  			
	list("a"=unlist(lapply(S,function(x) x$acb[1])),
		 "Theta"=unlist(lapply(S,function(x) .getAngle(x$angles[1]))),
		 "s"=unlist(lapply(S,function(x) x$acb[idx[1]]/x$acb[idx[2]])))
		 
}

#' Spatial histogram
#'
#' Characteristics for a spatial histogram
#'
#' Based on the estimated joint distribution after unfolding the function replicates
#' the entries of the vector \code{breaks}, see \code{\link{setbreaks}}, equal to the
#' number of estimated count data for each class. For an example please see the file 'unfold.R'.
#'
#' @param H 	  	trivariate (output) histogram, estimated by \code{\link{unfold}}
#' @param breaks 	breaks vector from \code{\link{setbreaks}}
#' 
#' @return list of the following entries: either major or minor semi-axis length \code{a},
#'         (polar) angle \code{Theta} and shape factor \code{s}, see \code{\link{parameters3d}}
#' 
#' @author M. Baaske
#' @rdname parameterEstimates
#' @export
parameterEstimates <- function(H,breaks) {
	list("a"=unlist(lapply(1:(length(breaks$size)-1),
							function(i) rep(breaks$size[i],sum(H[i,,])))),
		 "Theta"=unlist(lapply(1:(length(breaks$angle)-1),
							function(i) rep(breaks$angle[i],sum(H[,i,])))),
		 "s"=unlist(lapply(1:(length(breaks$shape)-1),
							function(i) rep(breaks$shape[i],sum(H[,,i])))))
}

#' Trivariate histogram
#'
#' 3D plot of a trivariate histogram
#' 
#' The (estimated spatial) joint size-shape-orientation distribution is plotted in a 3D
#' histogram with corresponding axes. The axes intersect in the first class number.
#' The ball volumes visualize the relative frequencies of count data for each class which
#' can be scaled by the user in order to make the spheres non-overlapping. Balls within the
#' same size class have the same color. For an example please see the file 'almmc.R'.
#'
#' @param A 		3D array of count data (estimated histogram), see \code{\link{unfold}}
#' @param main 		main title of the plot
#' @param scale 	scaling factor for non-overlapping spheres
#' @param col 		vector of color values repeatedly used for all classes
#' @param ... 		graphical parameters passed to \code{rgl::spheres3d}
#' 
#' @return 			NULL
#' 
#' @author M. Baaske
#' @rdname trivarHist
#' @export
trivarHist <- function(A, main = paste("Trivariate Histogram"), scale = 0.5, col, ...) {
  if (requireNamespace("rgl", quietly=TRUE)) {
	N <- sum(A)
	if(missing(col))
	 col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
	pos <- do.call(rbind,lapply(seq(1:dim(A)[1]),
					function(i) {
						X <- which(A[i,,]!=0,arr.ind=TRUE)
						cbind(X,rep(i,nrow(X)))
					}))
	## scaling of balls
	sz <- apply(pos, 1,function(x) scale*A[x[3],x[1],x[2]]/max(A))

	# coloring the balls
	xt <- as.vector(table(pos[,3]))
	cols2 <- rep(col,length.out=dim(A)[1])
	ncols <- lapply(1:dim(A)[1], function(i) rep(cols2[i], length.out=xt[i]))

	#plot spheres
	rgl::spheres3d(pos,radius=sz,col=unlist(ncols),...)
	rgl::axes3d(c('x','y','z'), pos=c(1,1,1), tick=FALSE)
	rgl::title3d(main,'',"orientation","shape","size")
 } else {
	 stop("Please install 'rgl' package from CRAN repositories before running this function.")
 }
}
