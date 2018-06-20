/**
 * simSphere.cpp
 *
 *  Created on: 07.05.2015
 *      Author: franke
 */

#include "SimSphere.h"

//static locals
static SEXP sphere_type_tag;
static int PL = 0;

#define MAX_ITER 100

SEXP convert_R_SphereSystem(STGM::Spheres& spheres, STGM::CBox3 &box);
SEXP convert_R_Circles(STGM::Intersectors<STGM::CSphere>::Type& objects);
STGM::CSphere convert_C_Sphere(SEXP R_sphere);
STGM::Spheres convert_C_Spheres(SEXP R_spheres);

R_Calldata buildRCallSpheres(SEXP R_param,SEXP R_cond) {
  int nprotect=0;
  R_Calldata d = Calloc(1,R_Calldata_s);

  d->call = R_NilValue;
  PROTECT(d->fname = getListElement( R_cond, "rdist")); ++nprotect;
  PROTECT(d->rho   = getListElement( R_cond, "rho"  )); ++nprotect;
  PROTECT(d->args  = getListElement( R_param,"radii")); ++nprotect;
  PROTECT(d->label = getListElement( R_cond, "label"));  ++nprotect;

  /* radii distribution */
  const char *ftype = CHAR(STRING_ELT(d->fname, 0));
  if ( !std::strcmp( ftype, "rlnorm") ||
       !std::strcmp( ftype, "rbeta" ) ||
       !std::strcmp( ftype, "rgamma") ||
       !std::strcmp( ftype, "runif" ) ||
       !std::strcmp( ftype, "const" ))
  {;
  } else {
    PROTECT(d->call = getCall(d->fname,d->args,d->rho)); ++nprotect;
  }
  d->nprotect = nprotect;
  d->isPerfect = asLogical(getListElement( R_cond, "perfect" ));
  return d;
}

void _free_spheres(STGM::CBoolSphereSystem *sp){
  if(!sp) return;
  sp->~CBoolSphereSystem();
  Free(sp);
}

void _sphere_finalizer(SEXP Rextp)
{
    checkPtr(Rextp, sphere_type_tag);
    STGM::CBoolSphereSystem *sp = (STGM::CBoolSphereSystem *)(R_ExternalPtrAddr(Rextp));
    _free_spheres(sp);
    R_ClearExternalPtr(Rextp);
}


SEXP FinalizeSphereSystem(SEXP ext) {
  _sphere_finalizer(ext);
  return R_NilValue;
}

SEXP CreateSpherePointer(STGM::CBoolSphereSystem *sp) {
  sphere_type_tag = install("SphereSystem_TAG");
  SEXP Rextp=R_NilValue;
  PROTECT(Rextp=R_MakeExternalPtr((void*) sp, sphere_type_tag, R_NilValue));
  R_RegisterCFinalizerEx(Rextp, (R_CFinalizer_t) _sphere_finalizer, TRUE);

  UNPROTECT(1);
  return Rextp;
}

STGM::CBoolSphereSystem * InitSphereSystem(SEXP R_param, SEXP R_cond) {
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  double lam   = asReal(AS_NUMERIC( getListElement( R_param, "lam")));

  /*  print level */
  PL = asInteger(getListElement( R_cond,"pl"));

  /* simulation box */
  STGM::CBox3 box(boxX,boxY,boxZ);

  /* set up sphere system */
  STGM::CBoolSphereSystem *sp = (STGM::CBoolSphereSystem*)Calloc(1,STGM::CBoolSphereSystem);

  try {
      new(sp)STGM::CBoolSphereSystem(box,lam);
  } catch(...) {
      error(_("InitSpheroidSystem(): Memory allocation error for sphere system."));
  }

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

  double m[3] = {m_box.m_size[0]+m_box.m_low[0],
                 m_box.m_size[1]+m_box.m_low[1],
                 m_box.m_size[2]+m_box.m_low[2]};

  /* loop over all */
  for (size_t niter=0;niter<num; niter++) {
      STGM::CVector3d center(runif(0.0,1.0)*m[0],runif(0.0,1.0)*m[1],runif(0.0,1.0)*m[2]);
      m_spheres.push_back( STGM::CSphere(center, f(), m_spheres.size()+1,label));
  }
}

void STGM::CBoolSphereSystem::simSpheresPerfect(double mx, double sdx, const char *label, int perfect) {
  int nTry=0, k=0;
  double p[4],sdx2=SQR(sdx), mu=0,r=0;;

  if(perfect) {
   cum_prob_k(mx,sdx2,m_box.m_up[0],m_box.m_up[1],m_box.m_up[2],p,&mu);
  } else mu = m_box.volume();

  /* get Poisson parameter */
  while(num==0 && nTry<MAX_ITER) {
        num = rpois(mu*m_lam);
        ++nTry;
  }
  m_spheres.reserve(num);

  if(PL>100) {
     if(perfect)
      Rprintf("Spheres (perfect) simulation, bivariate lognormal length/shape: \n");
      Rprintf("\t size distribution: %f %f %f \n",  mx,sdx,mu);
      Rprintf("\t cum sum of probabilities: %f, %f, %f, %f \n",p[0],p[1],p[2],p[3]);
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



void STGM::CBoolSphereSystem::simSphereSys(R_Calldata d)
{
   GetRNGstate();

   if(isNull(d->call)) {
       /* get arguments */
       double p1=REAL_ARG_LIST(d->args,0),p2=0;
       const char *fname = CHAR(STRING_ELT(d->fname,0));

       /* ACHTUNG: 'const' function braucht 2 Argumente */
       if(std::strcmp(fname, "const" ))
         p2=REAL_ARG_LIST(d->args,1);

       // set spheroid label
       const char *label = translateChar(asChar(d->label));

       if(!std::strcmp(fname, "rlnorm")) {
    	   simSpheresPerfect(p1,p2,label,d->isPerfect);
       } else {
           R_rndGen_t<rdist2_t> rrandom(p1,p2,fname);

           /* simulate with R's random generating functions */
           simSpheres<R_rndGen_t<rdist2_t> >(rrandom,label);
       }
   } else {
       /* eval R call for user defined radii distribution */
       const char *label = translateChar(asChar(d->label));
       R_eval_t<double> reval(d->call, d->rho);
       simSpheres<R_eval_t<double> &>(reval,label);
   }
   PutRNGstate();
}


SEXP GetSphereSystem(SEXP ext)
{
  checkPtr(ext, sphere_type_tag);

  STGM::CBoolSphereSystem *sp = (STGM::CBoolSphereSystem *)(getExternalPtr(ext));
  SEXP R_spheres = R_NilValue;
  STGM::Spheres &spheres = sp->refObjects();
  PROTECT(R_spheres = convert_R_SphereSystem(spheres,sp->box()));

  setAttrib(R_spheres, install("eptr"), ext);
  SET_CLASS_NAME(R_spheres,"spheres");

  UNPROTECT(1);
  return R_spheres;
}

SEXP SetupSphereSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond)
{
  SEXP R_var, R_ptr;
  PROTECT(R_var = getVar(R_vname,R_env));
  PROTECT(R_ptr = getAttrib(R_var,install("eptr")));

  if(isNull(R_ptr) || isNullPtr(R_ptr,sphere_type_tag)) {
     R_ptr = CreateSpherePointer( InitSphereSystem(R_param,R_cond) );
     if(PL>100) {
       Rprintf("setting pointer to %p \n",R_ExternalPtrAddr(R_ptr));
     }
  }

  STGM::CBoolSphereSystem *sp = (STGM::CBoolSphereSystem *)(getExternalPtr(R_ptr));
  sp->refObjects() = convert_C_Spheres(R_var);

  setAttrib(R_var, install("eptr"), R_ptr);
  UNPROTECT(2);
  return R_ptr;
}

SEXP SphereSystem(SEXP R_param, SEXP R_cond)
{
  STGM::CBoolSphereSystem *sp = InitSphereSystem(R_param,R_cond);
  R_Calldata cdata = buildRCallSpheres(R_param,R_cond);

  if(PL>100) Rprintf("Simulate... \n");
  sp->simSphereSys(cdata);

  if(PL>100) Rprintf("Simulated %d spheres.\n", sp->getNumSpheres());

  SEXP R_spheres=R_NilValue;
  if(PL>100) {
    Rprintf("Convert... \n");
    STGM::Spheres &spheres = sp->refObjects();
    PROTECT(R_spheres = convert_R_SphereSystem(spheres, sp->box()  ));
  } else {
    PROTECT(R_spheres = allocVector(VECSXP,0));  /* return empty list */
  }

  SEXP Rextp = R_NilValue;
  PROTECT(Rextp = CreateSpherePointer(sp));

  setAttrib(R_spheres, install("eptr"), Rextp);
  SET_CLASS_NAME(R_spheres,"sphere");

  deleteRCall(cdata);
  UNPROTECT(2);
  return R_spheres;
}


SEXP SimulateSpheresAndIntersect(SEXP R_param, SEXP R_cond, SEXP R_n) {
  STGM::CBoolSphereSystem *sp = InitSphereSystem(R_param,R_cond);
  R_Calldata cdata = buildRCallSpheres(R_param,R_cond);

  if(PL>100)
    Rprintf("Simulate and intersect... \n");

  sp->simSphereSys(cdata);
  deleteRCall(cdata);

  STGM::Intersectors<STGM::CSphere>::Type objects;
  int intern = 0;
  if(!isNull(getListElement(R_cond,"intern")))
   intern = asInteger(AS_INTEGER(getListElement(R_cond,"intern")));
  double dz = 0.0;
  if(!isNull(getListElement(R_cond,"dz")))
   dz = asReal(AS_NUMERIC(getListElement(R_cond,"dz")));
  else error(_("intersection is set to zero"));

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane(n,dz);
  //* intersection */
  sp->IntersectWithPlane(objects,plane,intern);

  //Rprintf("objects: %d , pl: %d, dz: %f \n",objects.size(),PL,dz);

  SEXP R_circles = R_NilValue;
  if(PL==10) {
    /* return radii only */
    PROTECT(R_circles = allocVector(REALSXP,objects.size()));
    for(size_t k=0;k<objects.size();k++)
      REAL(R_circles)[k] = objects[k].getCircle().r();
  } else {
    /* return full circle object */
    PROTECT(R_circles = convert_R_Circles(objects));
  }

  _free_spheres(sp);
  UNPROTECT(1);
  return R_circles;
}

SEXP IntersectSphereSystem(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_intern) {
  checkPtr(ext, sphere_type_tag);

  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane( n , asReal(R_z));

  STGM::CBoolSphereSystem *sp = static_cast<STGM::CBoolSphereSystem *>(getExternalPtr(ext));
  if(PL>100) Rprintf("Intersect with plane: %d \n", sp->refObjects().size());

  STGM::Intersectors<STGM::CSphere>::Type objects;
  int intern = asInteger(AS_INTEGER(R_intern));
  sp->IntersectWithPlane(objects,plane,intern);

  return convert_R_Circles(objects);
}


void STGM::CBoolSphereSystem::IntersectWithPlane(STGM::Intersectors<STGM::CSphere>::Type &objects, STGM::CPlane &plane, int intern)
{
  /// Intersect only objects fully inside the observation window
   int i=0,j=0;
   int l =plane.idx();
   switch(l) {
      case 0: i=1; j=2; break; // YZ
      case 1: i=0; j=2; break; // XZ
      case 2: i=0; j=1; break; // XY
   }

   // assume left-down corner is origin of box
   CWindow win(m_box.m_size[i],m_box.m_size[j]);

   for(size_t i=0; i<m_spheres.size(); ++i) {
       STGM::Intersector<STGM::CSphere> intersector( m_spheres[i], plane, m_box.m_size);
        if(intersector.FindIntersection()) {
          if(intersector.getCircle().isInWindow(win))
              objects.push_back( intersector );
        }
   }
}

STGM::CSphere convert_C_Sphere(SEXP R_sphere) {
  SEXP R_ctr;
  int interior=1;
  const char *label = "N";

  PROTECT(R_ctr = AS_NUMERIC( getListElement( R_sphere, "center")));

  int id = asInteger(AS_INTEGER( getListElement( R_sphere, "id")));
  double r = asReal(AS_NUMERIC(getListElement(R_sphere, "r")));

  if(!isNull(getAttrib(R_sphere, install("label"))))
    label = translateChar(asChar(getAttrib(R_sphere, install("label"))));

  if(!isNull(getAttrib(R_sphere, install("interior"))))
    interior = asLogical(getAttrib(R_sphere, install("interior")));

  UNPROTECT(1);
  return STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id,label,interior);
}


SEXP convert_R_SphereSystem(STGM::Spheres& spheres, STGM::CBox3 &box) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, spheres.size()) );

  STGM::CVector3d n(0,0,1);
  const STGM::LateralPlanes &planes = box.getLateralPlanes();

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<spheres.size();k++)
  {
    STGM::CSphere &sphere = spheres[k];

    STGM::Intersector<STGM::CSphere> intersector(sphere, box.m_size);
    Rboolean interior = (Rboolean) TRUE;
    for(size_t j=0; j<planes.size() ; ++j) {
         if( intersector(planes[j])) {
           interior = (Rboolean) FALSE;
           break;
         }
    }

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, 3));

    REAL(R_center)[0]=sphere.center()[0];
    REAL(R_center)[1]=sphere.center()[1];
    REAL(R_center)[2]=sphere.center()[2];

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

    PROTECT(R_names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(R_names, 0, mkChar("id"));
    SET_STRING_ELT(R_names, 1, mkChar("center"));
    SET_STRING_ELT(R_names, 2, mkChar("r"));

    setAttrib(R_tmp, R_NamesSymbol, R_names);
    setAttrib(R_tmp, install("label"), mkString(sphere.label()) );
    setAttrib(R_tmp, install("interior"), ScalarLogical(interior));
    setAttrib(R_tmp, install("area"), ScalarReal(sphere.projectionArea()));

    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(3);
  }

  UNPROTECT(1);
  return R_resultlist;
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres) {
  SEXP R_tmp, R_ctr;
  int id=0,
      N=length(R_spheres);
  STGM::Spheres spheres;
  spheres.reserve(N);

  double r=0;
  int interior=1;
  const char *label = "N";
  for(int i=0; i<N; i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheres,i));
      PROTECT(R_ctr = AS_NUMERIC( getListElement( R_tmp, "center")));
      id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
      r = asReal(getListElement( R_tmp, "r"));

      if(!isNull(getAttrib(R_tmp, install("label"))))
        label = translateChar(asChar(getAttrib(R_tmp, install("label"))));

      if(!isNull(getAttrib(R_tmp, install("interior"))))
       interior = asLogical(getAttrib(R_tmp, install("interior")));

      spheres.push_back(STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id,label,interior));
      UNPROTECT(2);
  }

  return spheres;
}


SEXP convert_R_Circles(STGM::Intersectors<STGM::CSphere>::Type & objects) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<objects.size();k++)
  {
     STGM::CCircle3 &circle = objects[k].getCircle();

     PROTECT(R_tmp = allocVector(VECSXP,ncomps));
     PROTECT(R_center = allocVector(REALSXP, 3));

     REAL(R_center)[0]=circle.center()[0];
     REAL(R_center)[1]=circle.center()[1];
     REAL(R_center)[2]=circle.center()[2];

     SET_VECTOR_ELT(R_tmp,0,ScalarInteger(circle.Id()));
     SET_VECTOR_ELT(R_tmp,1,R_center);
     SET_VECTOR_ELT(R_tmp,2,ScalarReal(circle.r()));

     PROTECT(R_names = allocVector(STRSXP, ncomps));
     SET_STRING_ELT(R_names, 0, mkChar("id"));
     SET_STRING_ELT(R_names, 1, mkChar("center"));
     SET_STRING_ELT(R_names, 2, mkChar("r"));

     setAttrib(R_tmp, R_NamesSymbol, R_names);
     SET_VECTOR_ELT(R_resultlist,k,R_tmp);
     UNPROTECT(3);
   }

   UNPROTECT(1);
   return R_resultlist;
 }
