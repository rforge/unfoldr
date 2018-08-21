\dontrun{

options(par.unfoldr=2L)
breaks <- setbreaks(c(6,5,6),maxSize=0.37,kap=1.25)
breaks

P <- coefficientMatrixSpheroids(breaks,check=FALSE)
c(min(P),max(P),sum(P))

}
