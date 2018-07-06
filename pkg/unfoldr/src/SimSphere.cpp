/**
 * simSphere.cpp
 *
 *  Created on: 07.05.2015
 *      Author: franke
 */

#include "SimSphere.h"

//static locals
static int PL = 0;

#define MAX_ITER 100

SEXP convert_R_SphereSystem(STGM::Spheres& spheres, STGM::CBox3 &box);
SEXP convert_R_Circles(STGM::Intersectors<STGM::CSphere>::Type& objects);
STGM::Spheres convert_C_Spheres(SEXP R_spheres);

void _free_spheres(STGM::CBoolSphereSystem *sp){
  if(!sp) return;
  sp->~CBoolSphereSystem();
  Free(sp);
}

STGM::CSphere convert_C_Sphere(SEXP R_sphere) {
  SEXP R_ctr;
  int interior=1;
  const char *label = "N";

  PROTECT(R_ctr = AS_NUMERIC( getListElement( R_sphere, "center")));

  int id = INTEGER(AS_INTEGER( getListElement( R_sphere, "id")))[0];
  double r = REAL(AS_NUMERIC(getListElement(R_sphere, "r")))[0];

  if(!isNull(getAttrib(R_sphere, install("label"))))
    label = translateChar(asChar(getAttrib(R_sphere, install("label"))));
  else {
  	label = "N";
  	warning(_("Setting `label` to `N` because it is undefined."));
  }

  if(!isNull(getAttrib(R_sphere, install("interior"))))
    interior = LOGICAL(getAttrib(R_sphere, install("interior")))[0];
  else {
	interior = 0;
	warning(_("Setting `interior` to zero because it is undefined."));
  }
  UNPROTECT(1);
  return STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id,label,interior);
}

STGM::CBoolSphereSystem * allocSphereSystem(STGM::CBox3 &box, double lam, int perfect) {
	/* set up sphere system */
	STGM::CBoolSphereSystem *sp = (STGM::CBoolSphereSystem*)Calloc(1,STGM::CBoolSphereSystem);
	try {
   	 new(sp)STGM::CBoolSphereSystem(box,lam,perfect);
	} catch(...) {
	 Rf_error(_("InitSpheredSystem(): Memory allocation error for sphere system."));
	}
	return sp;
}

STGM::CBoolSphereSystem * InitSphereSystem(SEXP R_cond) {
  SEXP R_box = R_NilValue;
  PROTECT(R_box  = getListElement( R_cond, "box"));

  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  /*  print level */
  PL = INTEGER(AS_INTEGER(getListElement( R_cond,"pl")))[0];
  int perfect = INTEGER_POINTER(getListElement( R_cond,"perfect"))[0];
  double lam  = NUMERIC_POINTER( getListElement( R_cond, "lam"))[0];

  /* simulation box */
  STGM::CBox3 box(boxX,boxY,boxZ);
  STGM::CBoolSphereSystem * sp = allocSphereSystem(box,lam,perfect);
  UNPROTECT(1);
  return sp;
}


template< typename  F>
void STGM::CBoolSphereSystem::simSpheres(F f, const char *label) {
  int nTry = 0;
  while(num==0 && nTry<MAX_ITER) {
     num = rpois(m_box.volume()*m_lam);
     ++nTry;
  }
  m_spheres.reserve(num);
  if(PL>10){
   Rprintf("box volume: %f, lam %f, number of spheres to simulate %d \n",m_box.volume(),m_lam,num);
  }
  double m[3] = {m_box.m_size[0]+m_box.m_low[0],
                 m_box.m_size[1]+m_box.m_low[1],
                 m_box.m_size[2]+m_box.m_low[2]};

  /* loop over all */
  for (size_t niter=0;niter<num; niter++)
  {
      STGM::CVector3d center(runif(0.0,1.0)*m[0],runif(0.0,1.0)*m[1],runif(0.0,1.0)*m[2]);
      m_spheres.push_back( STGM::CSphere(center, f(), m_spheres.size()+1,label));
  }
}

void STGM::CBoolSphereSystem::simSpheresPerfect(double mx, double sdx, const char *label, int perfect) {
  int nTry=0, k=0;
  double p[4],sdx2=SQR(sdx),mu=0,r=0;;

  if(perfect) {
   cum_prob_k(mx,sdx2,m_box.m_up[0],m_box.m_up[1],m_box.m_up[2],p,&mu);
  } else mu = m_box.volume();

  /* get Poisson parameter */
  while(num==0 && nTry<MAX_ITER) {
        num = rpois(mu*m_lam);
        ++nTry;
  }
  m_spheres.reserve(num);

  if(PL>10) {
     if(perfect){
       Rprintf("Spheres (perfect) simulation, lognormal radius. \n");
       Rprintf("\t size distribution: mx=%f, sdx=%f, mu=%f \n",  mx,sdx,mu);
       Rprintf("\t cum sum of probabilities: %f, %f, %f, %f \n",p[0],p[1],p[2],p[3]);
     }
  }

  /* loop over all */
  if(perfect) {
	  for (size_t niter=0;niter<num; niter++) {
	        sample_k(p,k);
	        r = rlnorm(mx+k*sdx2,sdx);
	        STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
	                               runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
	                               runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));

	        m_spheres.push_back( STGM::CSphere(center, r, m_spheres.size()+1,label));
	    }
  } else {
	  for (size_t niter=0;niter<num; niter++) {
	        r = rlnorm(mx,sdx);
	        STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
	                               runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
	                               runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

	        m_spheres.push_back( STGM::CSphere(center, r, m_spheres.size()+1,label));
	  }
  }

}


void STGM::CBoolSphereSystem::simSphereSys(SEXP R_args, SEXP R_cond) {
  SEXP R_fname, R_label;
  PROTECT(R_fname = getListElement( R_cond, "rdist"));
  PROTECT(R_label = getListElement( R_cond, "label"));

  /* radii distribution */
   const char *ftype = CHAR(STRING_ELT(R_fname, 0));
   const char *label = translateChar(asChar(R_label));

   GetRNGstate();
   if (!std::strcmp( ftype, "rlnorm") ||
	   !std::strcmp( ftype, "rbeta" ) ||
	   !std::strcmp( ftype, "rgamma") ||
	   !std::strcmp( ftype, "runif" ) ||
	   !std::strcmp( ftype, "const" ))
   {
	   double p1 = REAL_ARG_LIST(R_args,0);
       double p2 = LENGTH(R_args) > 1 ? REAL_ARG_LIST(R_args,1) : 0;

       if(PL>10) Rprintf("sphere simulation parameters: %f, %f ( %s ) \n",p1,p2, ftype);

	   if(!std::strcmp(ftype, "rlnorm")) {
		  int isPerfect = 0;
		  SEXP R_exact = R_NilValue;
		  PROTECT(R_exact = getListElement( R_cond, "perfect" ));
		  if(!isNull(R_exact) && !isLogical(R_exact))
			 isPerfect = LOGICAL(getListElement( R_cond, "perfect" ))[0];
		  simSpheresPerfect(p1,p2,label,isPerfect);
		  UNPROTECT(1);

	   } else {
		  R_rndGen_t<rdist2_t> rrandom(p1,p2,ftype);
		  simSpheres<R_rndGen_t<rdist2_t> >(rrandom,label);
	   }

   } else {
	  if(PL>10) Rprintf("simulate acc. to user-defined distribution ... \n");
	  SEXP R_call, R_rho;
	  PROTECT(R_rho = getListElement( R_cond, "rho" ));
	  PROTECT(R_call = getCall(R_fname,R_args,R_rho));
	  R_eval_t<double> reval(R_call,R_rho);
	  simSpheres<R_eval_t<double> &>(reval,label);
	  UNPROTECT(2);
   }

   PutRNGstate();
   if(PL>10) Rprintf("simulated %d spheres.\n", m_spheres.size());
   UNPROTECT(2);
   return;
}


SEXP SphereSystem(SEXP R_param, SEXP R_cond)
{
  STGM::CBoolSphereSystem *sp = InitSphereSystem(R_cond);
  if(PL>10)
	Rprintf("Simulate... \n");

  sp->simSphereSys(R_param,R_cond);
  STGM::Spheres &spheres = sp->refObjects();

  SEXP R_spheres = R_NilValue;
  if(PL==10) {
      /* return radii only */
	  size_t num = sp->getNumSpheres();
	  Rprintf("number of spheres: %d \n", num);

      PROTECT(R_spheres = allocVector(REALSXP,num));
      for(size_t k=0; k<num; k++)
        REAL(R_spheres)[k] = spheres[k].r();

  } else {
      /* return full circle object */
      PROTECT(R_spheres = convert_R_SphereSystem(spheres, sp->box()));
      SET_CLASS_NAME(R_spheres,"sphere");
  }

  _free_spheres(sp);
  UNPROTECT(1);
  return R_spheres;
}


SEXP SimulateSpheresAndIntersect(SEXP R_param, SEXP R_cond) {
  int nprotect = 0;
  /* init */
  STGM::CBoolSphereSystem *sp = InitSphereSystem(R_cond);
  if(PL>10){
  	Rprintf("simulate sphere system...\n");
  }
  sp->simSphereSys(R_param,R_cond);
  STGM::Intersectors<STGM::CSphere>::Type objects;

  int intern = 0;
  SEXP R_intern = R_NilValue;
  PROTECT(R_intern = AS_INTEGER(getListElement(R_cond,"intern"))); ++nprotect;
  if(!isNull(getListElement(R_cond,"intern")))
   intern = INTEGER(R_intern)[0];
  else { warning(_("Undefined argument `intern`. Set to zero.")); }

  double dz = 0.0;
  SEXP R_dz = R_NilValue;
  PROTECT(R_dz = AS_NUMERIC(getListElement(R_cond,"dz"))); ++nprotect;
  if(!isNull(R_dz))
   dz = REAL(R_dz)[0];
  else warning(_("Intersection coordinate is set to zero"));

  SEXP R_n = R_NilValue;
  PROTECT(R_n = AS_NUMERIC(getListElement(R_cond,"nsect"))); ++nprotect;
  if(isNull(R_n))
   error(_("Normal vector defining intersecting plane is `Null`."));

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane(n,dz);
  //* intersection */
  sp->IntersectWithPlane(objects,plane,intern);

  if(PL>10){
	  Rprintf("Plane normal to: [%f %f %f] \n", n[0],n[1],n[2]);
	  Rprintf("Number of intersections: %d at %f (intern=%d) \n",objects.size(),dz,intern);
  }

  SEXP R_circles = R_NilValue;
  if(PL==10) {
    /* return radii only */
    PROTECT(R_circles = allocVector(REALSXP,objects.size())); ++nprotect;
    for(size_t k=0;k<objects.size();k++)
      REAL(R_circles)[k] = objects[k].getCircle().r();

  } else {
    /* return full circle object */
    PROTECT(R_circles = convert_R_Circles(objects)); ++nprotect;
  }

  _free_spheres(sp);
  UNPROTECT(nprotect);
  return R_circles;
}

SEXP IntersectSphereSystem(SEXP R_var, SEXP R_n, SEXP R_z, SEXP R_intern, SEXP R_env, SEXP R_pl) {
  SEXP R_S;
  PROTECT(R_S = getVar(R_var,R_env));
  PL = INTEGER(AS_INTEGER(R_pl))[0];

  SEXP R_box;
  PROTECT(R_box = getAttrib(R_S, install("box")));
  if(isNull(R_box))
    error(_("Sphere system is missing a simulation box."));

  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  SEXP R_lam = R_NilValue;
  PROTECT(R_lam = getAttrib(R_S,install("lam")));
  if(isNull(R_lam))
 	error(_("Intensity parameter must be given as an attribute."));
  double lam = REAL(AS_NUMERIC(R_lam))[0];

  SEXP R_perfect = R_NilValue;
  PROTECT(R_perfect = getAttrib(R_S,install("perfect")));
  if(isNull(R_perfect))
  	  error(_("Whether simulation was exact or not must be given as an attribute."));
  int perfect = INTEGER(AS_INTEGER(R_perfect))[0];

  STGM::CBox3 box(boxX,boxY,boxZ);
  STGM::CBoolSphereSystem * sp = allocSphereSystem(box, lam, perfect);

  sp->refObjects() = convert_C_Spheres(R_S);

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane( n , REAL(AS_NUMERIC(R_z))[0]);

  if(PL>10)
	Rprintf("Intersect with plane: %d \n", sp->refObjects().size());

  STGM::Intersectors<STGM::CSphere>::Type objects;
  int intern = 0;
  if(!isNull(R_intern))
   intern = INTEGER(AS_INTEGER(R_intern))[0];
  sp->IntersectWithPlane(objects,plane,intern);
  if(PL>10) Rprintf("successful intersection of spheres. \n");

  SEXP R_circles = R_NilValue;
  if(PL == 10) {
      /* return radii only */
      PROTECT(R_circles = allocVector(REALSXP,objects.size()));
      for(size_t k=0;k<objects.size();k++)
        REAL(R_circles)[k] = objects[k].getCircle().r();
  } else {
      /* return full circle object */
      PROTECT(R_circles = convert_R_Circles(objects));
  }

  _free_spheres(sp);
  UNPROTECT(5);
  return R_circles;
}


void STGM::CBoolSphereSystem::IntersectWithPlane(STGM::Intersectors<STGM::CSphere>::Type &objects,
												 STGM::CPlane &plane, int intern)
{
  /// Intersect only objects fully inside the observation window
   int i=0,j=0, l =plane.idx();
   switch(l) {
      case 0: i=1; j=2; break; // YZ
      case 1: i=0; j=2; break; // XZ
      case 2: i=0; j=1; break; // XY
   }

   // assume left-down corner is origin of box
   CWindow win(m_box.m_size[i],m_box.m_size[j]);

   if(intern) {
	   for(size_t i=0; i<m_spheres.size(); ++i) {
		   STGM::Intersector<STGM::CSphere> intersector( m_spheres[i], plane, m_box.m_size);
			if(intersector.FindIntersection()) {
			  if(intersector.getCircle().isInWindow(win))
				  objects.push_back( intersector );
			}
	   }
   } else {
	   for(size_t i=0; i<m_spheres.size(); ++i) {
	   		STGM::Intersector<STGM::CSphere> intersector( m_spheres[i], plane, m_box.m_size);
	   		if(intersector.FindIntersection())
	   		   objects.push_back( intersector );
	    }
   }
}

SEXP convert_R_SphereSystem(STGM::Spheres& spheres, STGM::CBox3 &box) {
  if(PL>10) Rprintf("converting spheres ... \n");

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, spheres.size()) );

  STGM::CVector3d n(0,0,1);
  const STGM::LateralPlanes &planes = box.getLateralPlanes();

  SEXP R_tmp, R_center;
  const char *nms[] = {"id", "center", "r", ""};

  for(size_t k=0;k<spheres.size();k++) {
    STGM::CSphere &sphere = spheres[k];

    /* does the sphere touch the bounding box? (inner/outer particle) */
    STGM::Intersector<STGM::CSphere> intersector(sphere, box.m_size);
    Rboolean interior = (Rboolean) TRUE;
    for(size_t j=0; j<planes.size() ; ++j) {
         if( intersector(planes[j])) {
           interior = (Rboolean) FALSE;
           break;
         }
    }

    PROTECT(R_tmp = mkNamed(VECSXP, nms));
    PROTECT(R_center = allocVector(REALSXP, 3));

    REAL(R_center)[0]=sphere.center()[0];
    REAL(R_center)[1]=sphere.center()[1];
    REAL(R_center)[2]=sphere.center()[2];

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

    setAttrib(R_tmp, install("label"), mkString(sphere.label()) );
    setAttrib(R_tmp, install("interior"), ScalarLogical(interior));
    setAttrib(R_tmp, install("area"), ScalarReal(sphere.projectionArea()));

    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(2);
  }

  UNPROTECT(1);
  return R_resultlist;
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres) {
  if(PL>10) Rprintf("converting %d spheres ... \n", length(R_spheres));

  SEXP R_tmp, R_ctr;
  STGM::Spheres spheres;
  spheres.reserve(length(R_spheres));

  int interior = 1;
  const char *label = "N";
  for(int i=0; i<length(R_spheres); i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheres,i));
      PROTECT(R_ctr = AS_NUMERIC(getListElement( R_tmp, "center")));

      if(!isNull(getAttrib(R_tmp, install("label"))))
        label = translateChar(asChar(getAttrib(R_tmp, install("label"))));
      else {
    	label = "N";
    	warning(_("Undefined attribute `label`. Set to 'N'."));
      }

      if(!isNull(getAttrib(R_tmp, install("interior"))))
       interior = LOGICAL(getAttrib(R_tmp, install("interior")))[0];
      else {
    	interior = 0;
    	warning(_("Undefined attribute `interior`. Set to zero."));
      }

      spheres.push_back(STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],
    		  REAL(getListElement( R_tmp, "r"))[0],
			  INTEGER(AS_INTEGER(getListElement( R_tmp, "id")))[0],
			  label,interior));

      UNPROTECT(2);
  }
  return spheres;
}


SEXP convert_R_Circles(STGM::Intersectors<STGM::CSphere>::Type & objects) {
  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp, R_center;
  const char *nms[] = {"id", "center", "r", ""};

  for(size_t k=0;k<objects.size();k++)
  {
     STGM::CCircle3 &circle = objects[k].getCircle();
     PROTECT(R_tmp = mkNamed(VECSXP, nms));
     PROTECT(R_center = allocVector(REALSXP, 3));

     REAL(R_center)[0]=circle.center()[0];
     REAL(R_center)[1]=circle.center()[1];
     REAL(R_center)[2]=circle.center()[2];

     SET_VECTOR_ELT(R_tmp,0,ScalarInteger(circle.Id()));
     SET_VECTOR_ELT(R_tmp,1,R_center);
     SET_VECTOR_ELT(R_tmp,2,ScalarReal(circle.r()));

     SET_VECTOR_ELT(R_resultlist,k,R_tmp);
     UNPROTECT(2);
   }

   UNPROTECT(1);
   return R_resultlist;
 }
