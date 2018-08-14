\dontrun{

## Unfolding of the joint size-shape-orientation
## distribution of intersection ellipses of
## the data set `data15p`

library(unfoldr)
library(parallel)
options(par.unfoldr=detectCores())

# load ellipse intersection parameters
data(data15p)
# construct section profiles
AC <- data.matrix(data15p[c("A","C")])/1000 # unit: micro meter	

# for prolates: selecting the minor semi-axis lengths
# independent of nomenclature, which is always sp$A
sp <- sectionProfiles(AC,as.numeric(unlist(data15p["alpha"])))
summary(sp$alpha)

# set number of bins for each parameter
bin <- c(14,16,18)
# unfold the joint distribution
ret <- unfold(sp,bin,kap=1.25)

## optional: histogram plot
# library(rgl)
# trivarHist(ret$N_V,main="Trivariate Histogram (3d estimated)",scale=1.2)
	
## show marginal distributions
breaks <- ret$breaks
paramEst <- parameterEstimates(ret$N_V,ret$breaks)

#pdf("unfoldData.pdf",width = 8, height = 10)
op <- par(mfrow = c(3, 1))
# size
hist(paramEst$a,
		main=expression(paste("Minor semi-axis ",hat(c))),
		breaks=breaks$size,
		right=FALSE,freq=FALSE,col="gray",
		xlab=expression(hat(c)),ylim=c(0,25))

# Theta
hist(paramEst$Theta[paramEst$Theta<max(breaks$angle)],
		main=expression(paste("Polar angle ", theta)),
		breaks=breaks$angle,col="gray",right=FALSE,freq=FALSE,
		xlab=expression(theta),ylim=c(0,2.25))

# shape
hist(paramEst$s,main=expression(paste("Shape ",hat(s))),
		breaks=breaks$shape,xlim=c(0,1),ylim=c(0,8),
		right=FALSE,freq=FALSE,col="gray",xlab=expression(hat(s)))
par(op)

}
