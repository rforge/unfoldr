# Comment: Simulate Poisson sphere system and intersect 
# 
# Author: bahama
###############################################################################

library(rgl)
library(plotrix)
library(unfoldr)

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

#################################################################
## `beta` distribution for radii
#################################################################

lam <- 100
# log normal size, constant shape, 
theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5),
			  "shape"=list("s"=0.5),
			  "orientation"=list("kappa"=1))

# simulation bounding box
box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))

## simulate and return full spheres system
## intersect with XZ plane and return full list of intersection profiles

## section plane xy
head(simPoissonSystem)

S <- simPoissonSystem(theta,lam,size="rlnorm",box=box,type="prolate",
		intersect="original",mu=c(0,1,0),n=c(0,0,1),dz=0,perfect=TRUE,pl=101)

open3d()
spheroids3d(S[1:1500], FALSE, TRUE, box=box, col=col)
