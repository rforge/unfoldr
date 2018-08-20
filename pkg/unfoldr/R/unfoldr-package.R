###############################################################################
# Author:  M. Baaske
# Date:	   2018-08-20	
# File:    unfoldr-package.R: 
# 
# Comment: General description of the package
# 
###############################################################################

#' Stereological Unfolding for Spheroidal Particles
#'
#' Stereological unfolding as implemented in this package consists of the estimation of the joint size-shape-orientation
#' distribution of spheroidal shaped particles based on the same measured quantities of corresponding vertical
#' section profiles. A single trivariate discretized version of the (stereological) integral equation in the case of prolate
#' and oblate spheroids is solved numerically by a variant of the well-known Expectation Maximization (EM) algorithm. In addition,
#' routines for estimating the empirical diameter distribution of spheres from planar sections (better known as the Wicksell's
#' corpuscle problem) is also implemented. The package also provides functions for the simulation of Poisson germ-grain
#' processes with either spheroids, spherocylinders or spheres as grains including functions for planar and vertical sections
#' and digitization of section profiles.
#' 
#' @docType package
#' @name unfoldr-package
#' @importFrom stats na.omit
#' @useDynLib unfoldr, .registration = TRUE, .fixes = "C_"
#' 
#' @references  
#'  \enumerate{
#'     \item Bene\eqn{\check{\textrm{s}}},
#' 		  V. and Rataj, J. Stochastic Geometry: Selected Topics Kluwer Academic Publishers, Boston, 2004
#' 	   \item Ohser, J. and Schladitz, K. 3D images of materials structures Wiley-VCH, 2009
#' 	   \item Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
#'     \item C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms.
#' 					Springer, Berlin, 2002. Zbl 0990.86007
#' 	   \item M{\"u}ller, A., Weidner, A., and Biermann, H. (2015). Influence of reinforcement
#'            geometry on the very high-cycle fatigue behavior of aluminum-matrix-composites.
#' 			   Materials Science Forum, 825/826:150-157 
#'  } 	
#' 
NULL

#' Intersection ellipses parameters
#' 
#' The data set consists of section profile parameters (assumed to come from prolate spheroids) of fitted ellipses
#' based on measured section particles of an aluminium matrix composite [5] from metallographic analysis. The data
#' set can be used to reconstruct the trivariate spatial (prolate) spheroid distribution. 
#' 
#' @docType data
#' @keywords datasets
#' @name data15p
#' @usage data(data15p)
#' @format A matrix of columns: \code{A} (major semi-axis length), \code{C} (minor semi-axis length),
#' 	\code{S=C/A} (shape factor), \code{alpha} (polar angle in intersecting plane) and coordinates in the
#'  intersection plane \code{x,y}.
#' 
#' @author M. Baaske
NULL