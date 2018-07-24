# TODO: Add comment
# 
# Author: bahama
###############################################################################

library(rgl)
library(unfoldr)

spheres <- function(spheres, box=NULL, draw.axis=FALSE, draw.box=TRUE,...) {
	xyz <- apply(sapply(spheres, "[[", "center"),1,function(x) x)
	sizes <- unlist(lapply(spheres,function(x) x$r))
	spheres3d(xyz,radius=sizes,...)
	
	if(draw.axis) {
		axes3d(c('x','y','z'), pos=c(1,1,1), tick=FALSE)
		#title3d('','','x','y','z')		
	}
	
	## draw box
	if(draw.box & !is.null(box)) {
		x <- box$xrange[2]; y <- box$yrange[2]; z <- box$zrange[2]
		c3d.origin <- translate3d(scale3d(cube3d(col="darkgray", alpha=0.1),x/2,y/2,z/2),x/2,y/2,z/2)
		shade3d(c3d.origin)
		axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
				expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
	}
}

col <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF") 

## beta distribution for radii
lam <- 50
theta <- list("size"=list("shape1"=1,"shape2"=10))

box <- list("xrange"=c(0,5),"yrange"=c(0,5),"zrange"=c(0,5))




## return full spheres system
S <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",pl=101)


spheres(S$S,box,draw.box=TRUE,color=col)

## only intersections
head(simPoissonSystem)
sp <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",
		intersect="only", pl=10)


##**
r <- 0.5
X <- runif(100) * (5+2*r)+(-1)-r
max(X)
range(X)

box <- list("xrange"=c(-1,4),"yrange"=c(-1.5,3.5),"zrange"=c(0,2))
S <- simPoissonSystem(theta,lam,size="rbeta",box=box,type="sphere",
		intersect="original", pl=101)

notIn <- sapply(S,function(x) !attr(x,"interior"))


mx <- do.call(rbind,lapply(S[notIn],function(x) x$center))
apply(mx,2,min)

spheres(S[notIn],box,TRUE,draw.box=TRUE,color=col)
In <- sapply(S,function(x) attr(x,"interior"))
spheres(S[In],box,draw.box=TRUE,color="gray")


length(S)
summary(S)
theta$size[[1]]/(theta$size[[1]]+theta$size[[2]])			# mean
