\dontrun{

lam <- 10
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

# cylinders of constant length
theta <- list("size"=list(0.25),
			  "shape"=list(0.5),
			  "orientation"=list("kappa"=1))

S <- simCylinderSystem(theta,lam,size="const", shape="const",
						orientation="rbetaiso",box=box,pl=101)
				
# cylinders of constant length with
# beta distributed radius	
theta <- list("size"=list(0.35),
              "shape"=list("a"=1,"b"=5),
			  "orientation"=list("kappa"=1.5))				
				
S <- simCylinderSystem(theta,lam,size="const", shape="rbeta",
				orientation="rbetaiso",box=box,pl=101)
				
# bivariate length-shape distribution
# possibly correlated 
param <- list("mx"=-1.0,"my"=-2.5, "sdx"=0.15,"sdy"=0.2,"rho"=0.0,"kappa"=1.0)
theta <- list("size"=list("mx"=param$mx,"sdx"=param$sdx,
						  "my"=param$my,"sdy"=param$sdy,
						  "rho"=param$rho),
			  "orientation"=list("kappa"=param$kappa),
			  "shape"=list())
		
S <- simCylinderSystem(theta,lam,size="rbinorm",orientation="rbetaiso",box=box,pl=101)	

## show cylinder system
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#cylinders3d(S, box, col=cols)
 
}
