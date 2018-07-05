/**
 * @file SimCylinder.h
 *
 *  @date: 09.05.2016
 *  @author: M. Baaske
 */

//#include <vector>

#include "SimCylinder.h"
#include "directions.h"

static int PL = 0;

using namespace std;

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_i+DIM*_j] = (M)[_i][_j];         \
} while(0)

#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_i+DIM*_j];         \
} while(0)

STGM::CCylinder convert_C_Cylinder(SEXP R_cyl);
STGM::Cylinders convert_C_Cylinders(SEXP R_cyls);

SEXP convert_R_Cylinder( STGM::CCylinder &cyl, STGM::LateralPlanes &planes , STGM::CBox3 &box);
SEXP convert_R_Cylinders( STGM::Cylinders &cyl, STGM::CBox3 &box);
SEXP convert_R_CylinderIntersections(STGM::Intersectors<STGM::CCylinder>::Type &objects);

void _free_cylinders(STGM::CCylinderSystem *sp){
  if(!sp) return;
  sp->~CCylinderSystem();
  Free(sp);
}

SEXP UpdateIntersectionsCylinder(SEXP R_cylinder, SEXP R_cluster_ids, SEXP R_box) {
    int nProtected=0;
    SEXP R_ret, R_clust_id;

    SEXP R_BoxX, R_BoxY, R_BoxZ;
    PROTECT( R_BoxX = AS_NUMERIC( getListElement( R_box, "xrange" ) ) ); ++nProtected;
    PROTECT( R_BoxY = AS_NUMERIC( getListElement( R_box, "yrange" ) ) ); ++nProtected;
    PROTECT( R_BoxZ = AS_NUMERIC( getListElement( R_box, "zrange" ) ) ); ++nProtected;

    STGM::CBox3 box(REAL(R_BoxX)[1],REAL(R_BoxY)[1],REAL(R_BoxZ)[1]);
    const std::vector<STGM::CPlane> &planes = box.getPlanes();

    PROTECT(R_clust_id = getListElement(R_cluster_ids,"id")); ++nProtected;
    PROTECT(R_ret = allocVector(INTSXP,length(R_clust_id))); ++nProtected;

    for(int k=0;k<length(R_clust_id);k++) {
        STGM::CCylinder sp = convert_C_Cylinder(VECTOR_ELT(R_cylinder,k));
        STGM::Intersector<STGM::CCylinder> intersector(sp , box.m_size );
        for(size_t j=0; j<planes.size() ; ++j) {
            if( intersector(planes[j])) {
              sp.interior()=0;
              break;
            }
        }
        INTEGER(R_ret)[k]=sp.interior();
    }
    UNPROTECT(nProtected);
    return R_ret;
}

STGM::CCylinderSystem * allocCylinderSystem(STGM::CBox3 &box, double lam,
		STGM::CVector3d &maxis, int perfect)
{
	/* set up sphere system */
	STGM::CCylinderSystem *sp = (STGM::CCylinderSystem*)Calloc(1,STGM::CCylinderSystem);
	try {
		 new(sp)STGM::CCylinderSystem(box,lam,maxis,perfect);
	} catch(...) {
	 Rf_error(_("InitSpheroidSystem(): Memory allocation error for cylinder system."));
	}
	return sp;
}

STGM::CCylinderSystem * InitCylinderSystem(SEXP R_param, SEXP R_cond) {
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  // set print level
  PL = INTEGER(getListElement( R_cond,"pl"))[0];
  int perfect = INTEGER_POINTER(getListElement( R_cond,"perfect"))[0];
  double lam = NUMERIC_POINTER(getListElement( R_param, "lam"))[0];

  /* set up cylinder system */
  STGM::CBox3 box(boxX,boxY,boxZ);
  STGM::CVector3d maxis( NUMERIC_POINTER( getListElement( R_cond, "mu")) );
  STGM::CCylinderSystem *sp = allocCylinderSystem(box,lam,maxis,perfect);

  UNPROTECT(1);
  return sp;
}

void STGM::CCylinderSystem::simCylinderSystem(SEXP R_param, SEXP R_cond) {
	SEXP R_fname, R_args, R_label;
	PROTECT(R_fname = getListElement( R_cond, "rdist"));
	PROTECT(R_label = getListElement( R_cond, "label"));

	int isPerfect = LOGICAL(getListElement( R_cond, "perfect" ))[0];
	const char *label = translateChar(asChar(R_label));

	if(TYPEOF(R_fname) != VECSXP){     		 /* actually a function */
		SEXP R_call, R_rho;
		PROTECT(R_args = getListElement( R_param,"rmulti"));
		PROTECT(R_rho  = getListElement( R_cond, "rho" ));
		PROTECT(R_call = getCall(R_fname,R_args,R_rho));

		// simulate
		simJoint(R_call, R_rho, label);
		UNPROTECT(3);

	} else {
		/*
		 * fname, args are lists of [size, shape, orientation],
		 *  args are parameters for each
		 */
		PROTECT(R_args = allocVector(VECSXP,3));
	    SET_VECTOR_ELT(R_args,0, getListElement( R_param,"size"));
	    SET_VECTOR_ELT(R_args,1, getListElement( R_param,"shape"));
	    SET_VECTOR_ELT(R_args,2, getListElement( R_param,"orientation"));

	    // distribution types
	    const char *ftype_size  = GET_NAME(R_fname,0);
	    const char *ftype_shape = GET_NAME(R_fname,1);
	    const char *ftype_dir   = GET_NAME(R_fname,2);

	    STGM::CCylinder::direction_type dtype;
	    if (!std::strcmp( ftype_dir, "runifdir")) {
	    	dtype = STGM::CCylinder::UNIFORM_D;
	    } else if(!std::strcmp( ftype_dir, "rbetaiso" )) {
	    	dtype = STGM::CCylinder::BETAISOTROP_D;
	    } else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {
	    	dtype = STGM::CCylinder::MISES_D;
	    } else {
	       error(_("Direction distribution type is not supported."));
	    }

	    if( !std::strcmp( ftype_size, "rbinorm")) {

	    	simBivariate(R_args,dtype,label,isPerfect);

	    } else {

		    rdist2_t rshape;
		    if ( !std::strcmp(ftype_shape, "rbeta" )) {
		       rshape = &rbeta;
		    } else if(!std::strcmp(ftype_shape,"const")) {
		       rshape = &rconst;
		    } else if(!std::strcmp(ftype_shape,"rgamma")) {
		       rshape = &rgamma;
		    } else if(!std::strcmp(ftype_shape,"runif")) {
		       rshape = &runif;
		    } else {
		    	error(_("Unknown shape distribution type."));
		    }

		    rdist2_t rsize;
			if(!std::strcmp(ftype_size, "const" )) {
				rsize = &rconst;
			} else if (!std::strcmp(ftype_size, "rbeta" )) {
				rsize = &rbeta;
			} else if(!std::strcmp(ftype_size, "rlnorm")) {
				rsize = &rlnorm;
			} else if(!std::strcmp(ftype_size, "rgamma")) {
				rsize = &rgamma;
			} else if(!std::strcmp(ftype_size, "runif" )) {
				rsize = &runif;
			} else  {
			   error(_("Unknown size distribution type."));
			}

			simUnivar(R_args,rsize,rshape,dtype,label);

	    }
	    UNPROTECT(1);

	}

	UNPROTECT(2);
	return;
}


SEXP CylinderSystem(SEXP R_param, SEXP R_cond)
{
	/* init */
	STGM::CCylinderSystem *sp = InitCylinderSystem(R_param,R_cond);

	/* simulate  */
	sp->simCylinderSystem(R_param,R_cond);
	if(PL>100) Rprintf("Simulated %d spheroids: %p \n",sp->refObjects().size(), sp);

	SEXP R_cylinders = R_NilValue;
	if(PL>100) {
		STGM::Cylinders &cylinders = sp->refObjects();
		Rprintf("Convert... \n");
		PROTECT(R_cylinders = convert_R_Cylinders( cylinders, sp->box() ));
	} else {
		PROTECT(R_cylinders = allocVector(VECSXP,0));  /* return empty list */
	}
	classgets(R_cylinders, mkString("cylinder"));

	UNPROTECT(1);
	return R_cylinders;
}



SEXP SimulateCylindersAndIntersect(SEXP R_param, SEXP R_cond, SEXP R_n){
  int nprotect = 0;
  STGM::CCylinderSystem *sp = InitCylinderSystem(R_param,R_cond);

  /* simulate  */
  sp->simCylinderSystem(R_param,R_cond);
  if(PL>100) Rprintf("Simulated %d cylinders: %p \n",sp->refObjects().size(), sp);

  int intern = 0;
  SEXP R_intern = R_NilValue;
  PROTECT(R_intern = getListElement(R_cond,"intern")); ++nprotect;
  if(!isNull(R_intern))
	intern = INTEGER(AS_INTEGER(R_intern))[0];

  double dz = 0.0;
  SEXP R_dz = R_NilValue;
  PROTECT(R_dz = AS_NUMERIC(getListElement(R_cond,"dz"))); ++nprotect;
  if(!isNull(R_dz))
	 dz = REAL(R_dz)[0];
  else warning(_("Intersection coordinate is set to zero"));

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane(n,dz);

  STGM::Intersectors<STGM::CCylinder>::Type objects;
  sp->IntersectWithPlane(objects,plane,intern);

  if(PL>100) {
	  Rprintf("Plane normal:  [%f %f %f] \n", n[0],n[1],n[2]);
	  Rprintf("Total objects: %d \n", objects.size());
  }

  _free_cylinders(sp);
  UNPROTECT(nprotect);
  return convert_R_CylinderIntersections(objects);
}

void STGM::CCylinderSystem::simJoint(SEXP R_call, SEXP R_rho, const char *label) {
     GetRNGstate();

     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
       num = rpois(m_box.volume()*m_lam);
       ++nTry;
     }
     m_cylinders.reserve(num);

     // set box constants
     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double *v=0,h=0,radius=0,theta=0, phi=0,r=0;

     CVector3d u;
     SEXP Reval = R_NilValue;

     int err = 0;
     for (size_t niter=0; niter<num; niter++) {
         Reval = R_tryEval(R_call,R_rho,&err);
         if(!err) {
            h=2.0*REAL(getListElement(Reval,"h"))[0]; 		// simulate half length of cylinder !
            r=REAL(getListElement(Reval,"radius"))[0];
            theta=REAL(getListElement(Reval,"theta"))[0];
            phi=REAL(getListElement(Reval,"phi"))[0];

            v=REAL(getListElement(Reval,"u"));
            u[0]=v[0]; u[1]=v[1];  u[2]=v[2];

            STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
            m_cylinders.push_back(STGM::CCylinder(center,u,h,r,theta,phi,radius,niter+1,label) );
         } else
           error(_("simJoint(): `try` error in user defined distribution function."));
     }
     PutRNGstate();
}

void STGM::CCylinderSystem::simBivariate(SEXP R_args, STGM::CCylinder::direction_type dtype, const char *label, int perfect) {
    GetRNGstate();

    double mx=REAL(getListElement(VECTOR_ELT( R_args, 0),"mx"))[0];
	double my=REAL(getListElement(VECTOR_ELT( R_args, 0),"my"))[0];
	double sdx=REAL(getListElement(VECTOR_ELT( R_args, 0),"sdx"))[0];
	double sdy=REAL(getListElement(VECTOR_ELT( R_args, 0),"sdy"))[0];
	double rho=REAL(getListElement(VECTOR_ELT( R_args, 0),"rho"))[0];
	double kappa=REAL(getListElement(VECTOR_ELT(R_args,2),"kappa"))[0];

    /* get Poisson parameter */
	double p[4], mu = 0, sdx2 = SQR(sdx);

	if(perfect) {
		// cumulative probabilities
		cum_prob_k(mx,sdx2,m_box.m_up[0],m_box.m_up[1],m_box.m_up[2],p,&mu);
	} else {
		mu = m_box.volume();
	}

	if(PL>100) {
		 Rprintf("Spheroids (exact) simulation: \n");
		 Rprintf("\t size distribution: `rbinorm`  with %f %f %f %f %f\n",mx,my,sdx,sdy,rho);
		 Rprintf("\t directional distribution: %d  with %f \n", dtype, kappa);
		 Rprintf("\t cum sum of probabilities: %f, %f, %f, %f \n",p[0],p[1],p[2],p[3]);
		 Rprintf("\t set label: %s to character: \n",label);
		 Rprintf("\n\n");
	}

	int nTry=0;
	while(num==0 && nTry<MAX_ITER) {
		 num = rpois(mu*m_lam);
		 ++nTry;
	}
	m_cylinders.reserve(num);

    CVector3d u;
    double r=0,x=0,y=0,h=0,radius=0,s=1,phi=0,theta=0;

    if(perfect) {

      	  for (size_t niter=0; niter<num; niter++)
      	  {
      		  rbinorm_exact(p,mx,sdx,my,sdy,rho,x,y);
      		  s=1.0/(1.0+std::exp(-y));
			  h=std::exp(x);	  		/* overall length sampled including caps*/
			  radius=0.5*h*s;
			  h-=2.0*radius;  			/* then h is only height (excluding caps) */

			  if(m_maxR<r) m_maxR=r;	/* store maximum radius */
			  if(!R_FINITE(r))
				warning(_("simCylinderSysBivariat(): Some NA/NaN, +/-Inf produced"));

			   switch(dtype) {		    /* sample orientation */
				   case 0:
					 runidir(u.ptr(),theta,phi); break;
				   case 1:
					 if(kappa<1e-8) {
					   u = (runif(0.0,1.0)<0.5) ? m_mu : -m_mu;
					 } else {
					   rOhserSchladitz(u.ptr(),m_mu.ptr(), kappa, theta, phi);
					 }
					 break;
				   case 2:
					 if(kappa<1e-8)
					   runidir(u.ptr(),theta,phi);
					 else
					   rVonMisesFisher(u.ptr(), m_mu.ptr(), kappa, phi);
					 break;
			  }

			  /* sample positions conditionally of radii distribution */
			  STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
									 runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
									 runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));

			  m_cylinders.push_back(STGM::CCylinder(center,u,h,radius,theta,phi,r,niter+1,label) );
      	  }

    } else {

     	  for (size_t niter=0; niter<num; niter++)
     	  {
     		  rbinorm(mx,sdx,my,sdy,rho,x,y);
     		  s=1.0/(1.0+std::exp(-y));
			  h=std::exp(x);	  		/* overall length sampled including caps*/
			  radius=0.5*h*s;
			  h-=2.0*radius;  			/* then h is only height (excluding caps) */

			  switch(dtype) {		    /* sample orientation */
				   case 0:
					 runidir(u.ptr(),theta,phi); break;
				   case 1:
					 if(kappa<1e-8) {
					   u = (runif(0.0,1.0)<0.5) ? m_mu : -m_mu;
					 } else {
					   rOhserSchladitz(u.ptr(),m_mu.ptr(), kappa, theta, phi);
					 }
					 break;
				   case 2:
					 if(kappa<1e-8)
					   runidir(u.ptr(),theta,phi);
					 else
					   rVonMisesFisher(u.ptr(), m_mu.ptr(), kappa, phi);
					 break;
			  }

			  /* sample positions conditionally of radii distribution */
			  STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
									 runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
									 runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

			  m_cylinders.push_back(STGM::CCylinder(center,u,h,radius,theta,phi,r,niter+1,label) );
     	  }
    }

    PutRNGstate();

}

void STGM::CCylinderSystem::simUnivar(SEXP R_args, rdist2_t rsize, rdist2_t rshape,
			STGM::CCylinder::direction_type dtype, const char *label)
{
	 GetRNGstate();
	 int nTry=0;
	 while(num==0 && nTry<MAX_ITER) {
		num = rpois(m_box.volume()*m_lam);
		++nTry;
	 }
	 m_cylinders.reserve(num);

	 // size
	 double p2 = 0;
	 double p1 = REAL(VECTOR_ELT(VECTOR_ELT( R_args, 0),0))[0];
	 if(LENGTH(VECTOR_ELT( R_args, 0)) > 0)
	  p2 = REAL(VECTOR_ELT(VECTOR_ELT( R_args, 0),1))[0];

	 // shape
	 double s=1, s2=0;
	 double s1 = REAL(VECTOR_ELT(VECTOR_ELT(R_args,1),0))[0];
	 if(LENGTH(VECTOR_ELT( R_args, 1)) > 0)
		s2 = REAL(VECTOR_ELT(VECTOR_ELT( R_args, 1),1))[0];

	 // direction
	 double kappa = REAL(VECTOR_ELT(VECTOR_ELT(R_args,2),0))[0];

	 if(PL>100) {
	  Rprintf("Run spheroids  simulation... \n");
	  Rprintf("\t size distribution: %f %f \n", p1,p2);
	  Rprintf("\t shape parameters: %f %f\n", s1, s2);
	  Rprintf("\t directional distribution: %d  with %f \n", dtype, kappa);
	  Rprintf("\t set label: %s to character: \n",label);
	 }

	 /* loop over all */
     CVector3d u;
     double h=0, radius=0, theta=0, phi=0;
     for (size_t niter=0; niter<num; niter++)  {
         h = rsize(p1,p2);
         s = rshape(s1,s2);			/* could be constant = 1 */
         radius = 0.5*h*s;
         h-=2.0*radius;

         /* direction */
         switch(dtype) {
             case 0:
               runidir(u.ptr(),theta,phi); break;
             case 1:
               if(kappa<1e-8) {
                 u = (runif(0.0,1.0)<0.5) ? m_mu : -m_mu;
               } else {
                 rOhserSchladitz(u.ptr(),m_mu.ptr(), kappa, theta, phi);
               }
               break;
             case 2:
               if(kappa<1e-8)
                 runidir(u.ptr(),theta,phi);
               else
                 rVonMisesFisher(u.ptr(), m_mu.ptr(), kappa, phi);
               break;
        }

         /* sample positions conditionally of radii distribution */
         STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
                                runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
                                runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

         m_cylinders.push_back( STGM::CCylinder(center,u,h,radius,theta,phi,0,niter+1,label) );

     }
     PutRNGstate();
}

SEXP IntersectCylinderSystem(SEXP R_var, SEXP R_n, SEXP R_dz, SEXP R_intern, SEXP R_env, SEXP R_pl)
{
	  int nprotect = 0;
	  SEXP R_S = R_NilValue;
	  PROTECT(R_S = getVar(R_var,R_env)); ++nprotect;

	  int intern = 0;
	  PROTECT(R_intern = AS_INTEGER(R_intern)); ++nprotect;
	  if(!isNull(R_intern))
		intern = INTEGER(R_intern)[0];

	  double dz = 0.0;
	  PROTECT(R_dz = AS_NUMERIC(R_dz)); ++nprotect;
	  if(!isNull(R_dz))
		dz = REAL(R_dz)[0];
	  else warning(_("Intersection coordinate is set to zero"));

	  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
	  STGM::CPlane plane( n , dz);

	  SEXP R_box = R_NilValue;
	  PROTECT(R_box = getAttrib(R_S, install("box"))); ++nprotect;
	  if(isNull(R_box))
   		error(_("Cylinder system is missing a simulation box."));

	  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
	  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
	  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));
	  STGM::CBox3 box(boxX,boxY,boxZ);

	  SEXP R_mu = R_NilValue;
	  PROTECT(R_mu = getAttrib(R_S,install("R_mu"))); ++nprotect;
	  if(isNull(R_mu))
		error(_("Main orientation direction must not be 'Null'."));
	  STGM::CVector3d mu(REAL(R_mu)[0],REAL(R_mu)[1],REAL(R_mu)[2]);

	  SEXP R_lam = R_NilValue;
	  PROTECT(R_lam = getAttrib(R_S,install("lam"))); ++nprotect;
	  if(isNull(R_lam))
	  	error(_("Intensity parameter must be given as an attribute."));
	  double lam = REAL(AS_NUMERIC(R_lam))[0];

	  SEXP R_perfect = R_NilValue;
	  PROTECT(R_perfect = getAttrib(R_S,install("perfect")));
	  if(isNull(R_perfect))
	    error(_("Whether simulation was exact or not must be given as an attribute."));
	  int perfect = INTEGER(AS_INTEGER(R_perfect))[0];
	  /* alloc */
	  STGM::CCylinderSystem *sp = allocCylinderSystem(box,lam,mu,perfect);
	  // do conversions
	  sp->refObjects() = convert_C_Cylinders(R_S);
	  if(PL>100) {
	     Rprintf("Intersect with plane: %d , %p \n", sp->refObjects().size(), sp);
	  }

	  STGM::Intersectors<STGM::CCylinder>::Type objects;
	  sp->IntersectWithPlane(objects,plane,intern);

	  if(PL>100){
	   Rprintf("Plane normal:  [%f %f %f] \n", n[0],n[1],n[2]);
	   Rprintf("Total objects: %d \n", objects.size());
	  }

	  return convert_R_CylinderIntersections(objects);
}

SEXP convert_R_CylinderIntersections(STGM::Intersectors<STGM::CCylinder>::Type &objects)
{
  int  nLoopProtected=0, ncomps=16, ncompsCircle=4,inWindow=1;
  SEXP R_result, R_center, R_minor, R_major, R_ipt0, R_ipt1,
       R_mPoint0, R_mPoint1, R_height, R_ab, R_rcaps, R_psi,
	   R_obj = R_NilValue, R_names = R_NilValue;

  int type  = 0;
  PROTECT(R_result = allocVector(VECSXP, objects.size()) );
  //CWindow win(cylsys->getBox().m_size[0],cylsys->getBox().m_size[1]);

  for(size_t i=0; i<objects.size(); ++i)
  {
      type = (int) objects[i].getType();
      //inWindow = (int) objects[i].getEllipse().isInWindow(win);

      if(type == STGM::ELLIPSE ||
         type == STGM::ELLIPSE_ARC ||
         type == STGM::ELLIPSE_SEGMENT ||
         type == STGM::CIRCLE_CAPS) {

          PROTECT(R_obj = allocVector(VECSXP, ncomps) ); ++nLoopProtected;
          PROTECT(R_names = allocVector(STRSXP, ncomps));  ++nLoopProtected;

          PROTECT(R_center = allocVector(REALSXP, 3) );  ++nLoopProtected;
          PROTECT(R_minor = allocVector(REALSXP, 3) );   ++nLoopProtected;
          PROTECT(R_major = allocVector(REALSXP, 3) );   ++nLoopProtected;
          PROTECT(R_ipt0 = allocVector(REALSXP, 3) );    ++nLoopProtected;
          PROTECT(R_ipt1 = allocVector(REALSXP, 3) );    ++nLoopProtected;
          PROTECT(R_mPoint0 = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_mPoint1 = allocVector(REALSXP, 3) ); ++nLoopProtected;

          PROTECT(R_ab = allocVector(REALSXP, 2) );      ++nLoopProtected;
          PROTECT(R_psi = allocVector(REALSXP, 2) );     ++nLoopProtected;
          PROTECT(R_rcaps = allocVector(REALSXP, 2) );   ++nLoopProtected;
          PROTECT(R_height = allocVector(REALSXP, 2) );  ++nLoopProtected;

          STGM::CEllipse3 &ellipse = objects[i].getEllipse();
          REAL(R_center)[0]=ellipse.center()[0];
          REAL(R_center)[1]=ellipse.center()[1];
          REAL(R_center)[2]=ellipse.center()[2];

          REAL(R_major)[0]=ellipse.majorAxis()[0];
          REAL(R_major)[1]=ellipse.majorAxis()[1];
          REAL(R_major)[2]=ellipse.majorAxis()[2];

          REAL(R_minor)[0]=ellipse.minorAxis()[0];
          REAL(R_minor)[1]=ellipse.minorAxis()[1];
          REAL(R_minor)[2]=ellipse.minorAxis()[2];

          REAL(R_ipt0)[0]=objects[i].ipt0[0];
          REAL(R_ipt0)[1]=objects[i].ipt0[1];
          REAL(R_ipt0)[2]=objects[i].ipt0[2];

          REAL(R_ipt1)[0]=objects[i].ipt1[0];
          REAL(R_ipt1)[1]=objects[i].ipt1[1];
          REAL(R_ipt1)[2]=objects[i].ipt1[2];

          REAL(R_ab)[0] = ellipse.a();
          REAL(R_ab)[1] = ellipse.b();
          REAL(R_psi)[0]= ellipse.psi()[0];
          REAL(R_psi)[1]= ellipse.psi()[1];

          /* no conversion to [0,pi/2] because of plotting */
          double phi = objects[i].getCylinder().phi();

          /* convert angle */
          //if(phi > M_PI_2) {
		  // if(phi <= M_PI) phi = M_PI-phi;
		  // else if(phi < 1.5*M_PI) phi = std::fmod(phi,M_PI);
		  // else phi = 2.0*M_PI-phi;
          //}

          REAL(R_mPoint0)[0]=objects[i].getCircle1().center()[0];
          REAL(R_mPoint0)[1]=objects[i].getCircle1().center()[1];
          REAL(R_mPoint0)[2]=objects[i].getCircle1().center()[2];

          REAL(R_mPoint1)[0]=objects[i].getCircle2().center()[0];
          REAL(R_mPoint1)[1]=objects[i].getCircle2().center()[1];
          REAL(R_mPoint1)[2]=objects[i].getCircle2().center()[2];

          REAL(R_rcaps)[0]= objects[i].getCircle1().r();
          REAL(R_rcaps)[1]= objects[i].getCircle2().r();

          SET_VECTOR_ELT(R_obj,2, R_center );
          SET_VECTOR_ELT(R_obj,3, R_major  );
          SET_VECTOR_ELT(R_obj,4, R_minor  );
          SET_VECTOR_ELT(R_obj,5, R_ipt0   );
          SET_VECTOR_ELT(R_obj,6, R_ipt1   );
          SET_VECTOR_ELT(R_obj,7, R_mPoint0);
          SET_VECTOR_ELT(R_obj,8, R_mPoint1);
          SET_VECTOR_ELT(R_obj,9, R_ab     );
          SET_VECTOR_ELT(R_obj,10,ScalarReal(phi));
          SET_VECTOR_ELT(R_obj,11,ScalarReal(ellipse.b()/ellipse.a()));
          SET_VECTOR_ELT(R_obj,12,R_psi    );
          SET_VECTOR_ELT(R_obj,13,R_rcaps  );
          SET_VECTOR_ELT(R_obj,14,ScalarInteger(objects[i].getSide()));
          SET_VECTOR_ELT(R_obj,15,ScalarInteger(inWindow));

          SET_STRING_ELT(R_names, 2, mkChar("center"));
          SET_STRING_ELT(R_names, 3, mkChar("major"));
          SET_STRING_ELT(R_names, 4, mkChar("minor"));
          SET_STRING_ELT(R_names, 5, mkChar("ipt0"));
          SET_STRING_ELT(R_names, 6, mkChar("ipt1"));
          SET_STRING_ELT(R_names, 7, mkChar("mPoint0"));
          SET_STRING_ELT(R_names, 8, mkChar("mPoint1"));
          SET_STRING_ELT(R_names, 9, mkChar("ab"));
          SET_STRING_ELT(R_names, 10, mkChar("phi"));
          SET_STRING_ELT(R_names, 11, mkChar("shape"));
          SET_STRING_ELT(R_names, 12, mkChar("psi"));
          SET_STRING_ELT(R_names, 13, mkChar("rcaps"));
          SET_STRING_ELT(R_names, 14, mkChar("pS"));
          SET_STRING_ELT(R_names, 15, mkChar("inW"));

      } else if(type ==STGM::CIRCLE ){
          PROTECT(R_obj = allocVector(VECSXP, ncompsCircle) ); ++nLoopProtected;
          PROTECT(R_names = allocVector(STRSXP, ncompsCircle));  ++nLoopProtected;
          PROTECT(R_mPoint0 = allocVector(REALSXP, 3) );       ++nLoopProtected;

          /* circle1 stores the circle as intersection of cylidner */
          REAL(R_mPoint0)[0]=objects[i].getCircle1().center()[0];
          REAL(R_mPoint0)[1]=objects[i].getCircle1().center()[1];
          REAL(R_mPoint0)[2]=objects[i].getCircle1().center()[2];

          SET_VECTOR_ELT(R_obj,2, R_mPoint0 );
          SET_VECTOR_ELT(R_obj,3,ScalarReal( objects[i].getCircle1().r() ));

          SET_STRING_ELT(R_names, 2, mkChar("mPoint0"));
          SET_STRING_ELT(R_names, 3, mkChar("radius"));

      }

      SET_VECTOR_ELT(R_obj,0, ScalarInteger( (int) objects[i].getCylinder().Id()));
      SET_VECTOR_ELT(R_obj,1, ScalarInteger( type ));

      SET_STRING_ELT(R_names, 0, mkChar("id"));
      SET_STRING_ELT(R_names, 1, mkChar("type"));

      setAttrib(R_obj, R_NamesSymbol, R_names);
      SET_VECTOR_ELT(R_result,i,R_obj);

      UNPROTECT(nLoopProtected);
      nLoopProtected=0;
  }

  UNPROTECT(1);
  return R_result;
}


void STGM::CCylinderSystem::IntersectWithPlane(STGM::Intersectors<STGM::CCylinder>::Type &objects,CPlane &plane, int intern)
{
  int i=0,j=0;
  switch(plane.idx()) {
        case 0: i=1; j=2; break; // YZ
        case 1: i=0; j=2; break; // XZ
        case 2: i=0; j=1; break; // XY
  }
  CWindow win(m_box.m_size[i],m_box.m_size[j]);

  for(size_t i=0; i<m_cylinders.size(); ++i)
  {
         STGM::Intersector<STGM::CCylinder> intersector(m_cylinders[i], plane, m_box.m_size);
         if(intersector.TestIntersection()) {
             intersector.FindIntersection();
             if(intersector.getType() == CIRCLE_CAPS ||
                intersector.getType() == CIRCLE) {
               if(intern) {
                 if((intersector.getCircle1()).isInWindow(win))
                   objects.push_back(intersector);
               } else objects.push_back(intersector);
             } else if(intersector.getType() == ELLIPSE ||
                       intersector.getType() == ELLIPSE_ARC ||
                       intersector.getType() == ELLIPSE_SEGMENT) {

                 if(intern) {
                   if((intersector.getEllipse()).isInWindow(win))
                     objects.push_back(intersector);
                 } else objects.push_back(intersector);
             }
         }

     }
}

STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  double radius = 0;
  const char *label = "N";

  if(!isNull(getAttrib(R_cyl, install("label"))))
    label = translateChar(asChar(getAttrib(R_cyl, install("label"))));
  if(!isNull(getAttrib(R_cyl, install("interior"))))
    interior = LOGICAL(getAttrib(R_cyl, install("interior")))[0];
  if(!isNull(getAttrib(R_cyl, install("radius"))))
    radius = REAL(getAttrib(R_cyl, install("radius")))[0];

  if(!std::strcmp(label,"F")) {
      SEXP R_ctr;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      STGM::CVector3d ctr(REAL(R_ctr)),u(0,0,1);

      UNPROTECT(1);
      return STGM::CCylinder(ctr,u,0,REAL(getListElement(R_cyl, "r"))[0],0,0,
                 radius, INTEGER(getListElement(R_cyl, "id"))[0], label, interior);

  } else {
      SEXP R_ctr, R_u, R_angles;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_cyl, "u")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_cyl, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)),u(REAL(R_u));
      UNPROTECT(3);

      return STGM::CCylinder(ctr,u, REAL(getListElement(R_cyl, "length"))[0],
                REAL(getListElement(R_cyl, "r"))[0], REAL(R_angles)[0],REAL(R_angles)[1],
                  radius, INTEGER(getListElement(R_cyl, "id"))[0], label, interior);
  }

}

STGM::Cylinders convert_C_Cylinders(SEXP R_cyls)
{
  SEXP R_cyl, R_ctr, R_u, R_angles;

  SEXP R_label = R_NilValue;
  PROTECT(R_label = getAttrib(R_cyls, install("label")));
  const char *label = "N";
  if(!isNull(R_label))
    label = translateChar(asChar(R_label));

  int interior = 1;
  double radius = 0;
  STGM::Cylinders cylinders;

  for(int i=0; i<length(R_cyls); i++) {
      PROTECT(R_cyl = VECTOR_ELT(R_cyls,i));

      if(!isNull(getAttrib(R_cyl, install("label"))))
        label = translateChar(asChar(getAttrib(R_cyl, install("label"))));
      else {
    	 label = "N";
    	 warning(_("Attribute `label` is Null. Set to `N`."));
      }
      if(!isNull(getAttrib(R_cyl, install("interior"))))
        interior = LOGICAL(getAttrib(R_cyl, install("interior")))[0];
      else {
    	  interior = 0;
    	  warning(_("Attribute `interior` is Null. Set to `0`. "));
      }
      if(!isNull(getAttrib(R_cyl, install("radius"))))
        radius = REAL(getAttrib(R_cyl, install("radius")))[0];
      else {
    	  radius = 0;
    	  warning(_("Attribute `radius` is Null. Set to `0`."));
      }

      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_cyl, "u")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_cyl, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
      STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

      cylinders.push_back(
    	STGM::CCylinder(ctr,u,
    		  REAL(getListElement(R_cyl, "length"))[0],
    		  REAL(getListElement(R_cyl, "r"))[0],
    		  REAL(R_angles)[0],REAL(R_angles)[1],
			  radius,
			  INTEGER(getListElement(R_cyl, "id"))[0],
			  label, interior));

      UNPROTECT(4);
  }

  UNPROTECT(1);
  return cylinders;
}


SEXP convert_R_Cylinder( STGM::CCylinder &cyl, STGM::LateralPlanes &planes, STGM::CBox3 &box) {
  int ncomps=9, dim=3;
  SEXP R_names, R_Cyl = R_NilValue;
  SEXP R_u, R_center,  R_origin0, R_origin1, R_angles, R_rotM;

  PROTECT(R_Cyl = allocVector(VECSXP,ncomps));
  PROTECT(R_u = allocVector(REALSXP, 3));
  PROTECT(R_center = allocVector(REALSXP, 3));
  PROTECT(R_origin0 = allocVector(REALSXP, 3));
  PROTECT(R_origin1 = allocVector(REALSXP, 3));
  PROTECT(R_angles = allocVector(REALSXP,2));
  PROTECT(R_rotM = allocMatrix(REALSXP,3,3));

  REAL(R_u)[0]=cyl.u()[0];
  REAL(R_u)[1]=cyl.u()[1];
  REAL(R_u)[2]=cyl.u()[2];

  REAL(R_center)[0]=cyl.center()[0];
  REAL(R_center)[1]=cyl.center()[1];
  REAL(R_center)[2]=cyl.center()[2];

  REAL(R_origin0)[0] = cyl.origin0()[0];
  REAL(R_origin0)[1] = cyl.origin0()[1];
  REAL(R_origin0)[2] = cyl.origin0()[2];

  REAL(R_origin1)[0] = cyl.origin1()[0];
  REAL(R_origin1)[1] = cyl.origin1()[1];
  REAL(R_origin1)[2] = cyl.origin1()[2];

  REAL(R_angles)[0]=cyl.theta();
  REAL(R_angles)[1]=cyl.phi();

  //for (size_t i = 0; i < 3; i++)
  //  for (size_t j = 0; j < 3; j++)
  //    REAL(R_rotM)[i + 3 *j] = (cyl.rotationMatrix())[i][j];

  const STGM::CMatrix3d &M = cyl.rotationMatrix();
  COPY_C2R_MATRIX(M,R_rotM,dim);

  // check intersection
  STGM::Intersector<STGM::CCylinder> intersector(cyl, box.m_size );
  Rboolean interior = (Rboolean) TRUE;
  for(size_t j=0; j<planes.size() ; ++j) {
       if( intersector(planes[j])) {
         interior = (Rboolean) FALSE;
         break;
       }
  }
  STGM::PointVector2d p;
  double area = cyl.delamProjection(p,20);

  SET_VECTOR_ELT(R_Cyl,0,ScalarInteger( (int) cyl.Id() ));
  SET_VECTOR_ELT(R_Cyl,1,R_center);
  SET_VECTOR_ELT(R_Cyl,2,R_origin0);
  SET_VECTOR_ELT(R_Cyl,3,R_origin1);
  SET_VECTOR_ELT(R_Cyl,4,ScalarReal(cyl.h())) ;
  SET_VECTOR_ELT(R_Cyl,5,R_u);
  SET_VECTOR_ELT(R_Cyl,6,ScalarReal(cyl.r()));
  SET_VECTOR_ELT(R_Cyl,7,R_angles);
  SET_VECTOR_ELT(R_Cyl,8,R_rotM);

  PROTECT(R_names = allocVector(STRSXP, ncomps));
  SET_STRING_ELT(R_names, 0, mkChar("id"));
  SET_STRING_ELT(R_names, 1, mkChar("center"));
  SET_STRING_ELT(R_names, 2, mkChar("origin0"));
  SET_STRING_ELT(R_names, 3, mkChar("origin1"));
  SET_STRING_ELT(R_names, 4, mkChar("length"));
  SET_STRING_ELT(R_names, 5, mkChar("u"));
  SET_STRING_ELT(R_names, 6, mkChar("r"));
  SET_STRING_ELT(R_names, 7, mkChar("angles"));
  SET_STRING_ELT(R_names, 8, mkChar("rotM"));

  setAttrib(R_Cyl, R_NamesSymbol, R_names);
  setAttrib(R_Cyl, install("radius"), ScalarReal(cyl.radius()));
  setAttrib(R_Cyl, install("label"), mkString(cyl.label()) );
  setAttrib(R_Cyl, install("interior"), ScalarLogical(interior));
  // here always 'delam' projection area
  setAttrib(R_Cyl, install("area"), ScalarReal(area));

  UNPROTECT(8);
  return R_Cyl;
}

SEXP convert_R_Cylinders( STGM::Cylinders &cyls, STGM::CBox3 &box) {
  Rprintf("Convert...%d  cylinders.\n", cyls.size());

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, cyls.size()) );
  // get lateral bounding planes
  STGM::LateralPlanes &planes = box.getLateralPlanes();

  for(size_t k=0;k<cyls.size();++k)
    SET_VECTOR_ELT(R_resultlist,k,convert_R_Cylinder(cyls[k],planes,box));

  UNPROTECT(1);
  return(R_resultlist);
}


SEXP CDigitizeCylinderIntersections(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_delta)
{

/*
  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane( n , asReal(R_z));

  //STGM::IntersectorCylinders objects;
  STGM::Intersectors<STGM::CCylinder>::Type objects;
  STGM::CCylinderSystem *cyls = static_cast<STGM::CCylinderSystem *>(getExternalPtr(ext));
  cyls->IntersectWithPlane(objects,plane,0);

  int nPix = (int) cyls->box().m_size[0]/asReal(R_delta);

  SEXP R_W = R_NilValue;
  PROTECT(R_W = allocMatrix(INTSXP,nPix,nPix));

  Rprintf("Digitize: nPix: %d, delta: %f \n",nPix,REAL(R_delta)[0]);
  STGM::digitize<STGM::CCylinder>(objects,INTEGER(R_W),nPix,asReal(R_delta));

  UNPROTECT(1);
  return R_W;
*/

  return R_NilValue;
}
