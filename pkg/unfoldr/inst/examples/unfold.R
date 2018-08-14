\dontrun{

## Comment: Trivariate unfolding of spheroid distribution
# library(unfoldr)

## set number of cpu cores (optional)
# library(parallel)
# options(par.unfoldr=detectCores())

## Intensity: mean number of spheroids per unit volume
lam <- 2500

## simulation parameters
theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5),
		      "shape"=list(0.5),"orientation"=list("kappa"=2))
## simualtion
set.seed(1234)

S <- simPoissonSystem(theta,lam,size="rlnorm",
		orientation="rbetaiso",box=list(c(0,5)),type="prolate",pl=101)

## unfolding
sp <- verticalSection(S,2.5)
ret <- unfold(sp,c(15,12,11),kap=1.25)
cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")

## plot 3d trivariate histogram of joint distribution
#trivarHist(ret$N_V,scale=0.9)

## aggregate data for hists
param3d <- parameters3d(S)
paramEst <- parameterEstimates(ret$N_V,ret$breaks)

## Marginal histograms of
## size (minor semi-axis), shape and orientation

# pdf("spheroidHist.pdf",width = 8, height = 10)
op <- par(mfrow = c(3, 2))
hist(param3d$a[param3d$a<max(ret$breaks$size)],
 main=expression(paste("3D Histogram ", c)),
 breaks=ret$breaks$size,col="gray",right=FALSE,freq=FALSE,xlab="c",ylim=c(0,25))

hist(paramEst$a,
 main=expression(paste("Estimated histogram ",hat(c))),
 breaks=ret$breaks$size,
 right=FALSE,freq=FALSE,col="gray",
 xlab=expression(hat(c)),ylim=c(0,25))

hist(param3d$Theta[param3d$Theta<max(ret$breaks$angle)],
 main=expression(paste("3D Histogram ", theta)),
 breaks=ret$breaks$angle,col="gray",right=FALSE,freq=FALSE,
 xlab=expression(theta),ylim=c(0,2))

hist(paramEst$Theta,
 main=expression(paste("Estimated Histogram ", hat(theta))),
 breaks=ret$breaks$angle,
 right=FALSE,freq=FALSE,col="gray",
 xlab=expression(hat(theta)),ylim=c(0,2))

hist(param3d$s,main=expression(paste("3D Histogram ", s)),
 col="gray",breaks=ret$breaks$shape,
 right=FALSE,freq=FALSE,xlab="s",ylim=c(0,10))

hist(paramEst$s,main=expression(paste("Estimated Histogram ", hat(s))),
 breaks=ret$breaks$shape,
right=FALSE,freq=FALSE,col="gray",xlab=expression(hat(s)),ylim=c(0,10))
par(op)
#dev.off()
	
}
