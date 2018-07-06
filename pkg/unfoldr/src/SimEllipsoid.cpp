/**
 *  @file SimEllipsoid.cpp
 *
 *  @date: 10.01.2014
 *  @author: M. Baaske
 */

#define MAX_ITER 100

#include "SimEllipsoid.h"
#include "directions.h"

//static locals
static int PL = 0;

using namespace std;

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_j+DIM*_i] = (M)[_i][_j];         \
} while(0)

#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_j+DIM*_i];         \
} while(0)

#define GET_OBJECT_CLASS(RS) translateChar(asChar(getAttrib( (RS), R_ClassSymbol)))

SEXP convert_R_EllipsoidSystem( STGM::CSpheroidSystem *S);
SEXP convert_R_Ellipses_all(STGM::Intersectors<STGM::CSpheroid>::Type &objects);
SEXP convert_R_Ellipses_trunc(STGM::Intersectors<STGM::CSpheroid>::Type &objects);

extern STGM::CSphere convert_C_Sphere(SEXP R_sphere);
extern STGM::CCylinder convert_C_Cylinder(SEXP R_cyl);
extern STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid);

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids);

void _free_spheroids(STGM::CSpheroidSystem *sp){
  if(!sp) return;
  sp->~CSpheroidSystem();
  Free(sp);
}

STGM::CSpheroidSystem * allocSpheroidSystem(STGM::CBox3 &box, double lam,
		STGM::CVector3d &maxis, STGM::CSpheroid::spheroid_type stype, int perfect)
{
	/* set up sphere system */
	STGM::CSpheroidSystem *sp = (STGM::CSpheroidSystem*)Calloc(1,STGM::CSpheroidSystem);
	try {
		 new(sp)STGM::CSpheroidSystem(box,lam,maxis,stype,perfect);
	} catch(...) {
	 Rf_error(_("InitSpheroidSystem(): Memory allocation error for spheroid system."));
	}
	return sp;
}


STGM::CSpheroidSystem * InitSpheroidSystem(SEXP R_cond) {
  /* get the box */
  SEXP R_box;
  PROTECT( R_box  = getListElement( R_cond, "box"));
  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

  SEXP R_spheroidType = R_NilValue;
  PROTECT( R_spheroidType = getListElement( R_cond, "stype" ) );
  const char* stype_str = CHAR( STRING_ELT( R_spheroidType, 0 ));

  STGM::CSpheroid::spheroid_type stype = STGM::CSpheroid::PROLATE;
  if( !std::strcmp("oblate",stype_str))
    stype = STGM::CSpheroid::OBLATE;

  // set print level
  PL = INTEGER(AS_INTEGER(getListElement( R_cond,"pl")))[0];
  int perfect = INTEGER_POINTER(getListElement( R_cond,"perfect"))[0];
  double lam = NUMERIC_POINTER(getListElement( R_cond, "lam"))[0];

  /* set up spheroid system */
  STGM::CBox3 box(boxX,boxY,boxZ);

  double *mu = NUMERIC_POINTER( getListElement( R_cond, "mu"));
  STGM::CVector3d maxis(mu[0],mu[1],mu[2]);
  STGM::CSpheroidSystem *sp = allocSpheroidSystem(box,lam,maxis,stype,perfect);

  UNPROTECT(2);
  return sp;
}


void STGM::CSpheroidSystem::simSpheroidSystem(SEXP R_args, SEXP R_cond) {
	SEXP R_fname, R_label;
	PROTECT(R_fname = getListElement( R_cond, "rdist"));
	PROTECT(R_label = getListElement( R_cond, "label"));

	int isPerfect = LOGICAL(getListElement( R_cond, "perfect" ))[0];
	const char *label = translateChar(asChar(R_label));

	GetRNGstate();
	if(TYPEOF(R_fname) != VECSXP){
		SEXP R_call, R_rho;
		PROTECT(R_rho  = getListElement( R_cond, "rho" ));
		PROTECT(R_call = getCall(R_fname,R_args,R_rho));

		// simulate
		simJoint(R_call, R_rho, label);
		UNPROTECT(2);

	} else {

	    // distribution types
	    const char *ftype_size  = GET_NAME(R_fname,0);
	    const char *ftype_shape = GET_NAME(R_fname,1);
	    const char *ftype_dir   = GET_NAME(R_fname,2);

	    STGM::CSpheroid::direction_type dtype;
	    if (!std::strcmp( ftype_dir, "runifdir")) {
	    	dtype = STGM::CSpheroid::UNIFORM_D;
	    } else if(!std::strcmp( ftype_dir, "rbetaiso" )) {
	    	dtype = STGM::CSpheroid::BETAISOTROP_D;
	    } else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {
	    	dtype = STGM::CSpheroid::MISES_D;
	    } else {
	       error(_("Direction distribution type is not supported."));
	    }

	    rdist2_t rshape;
	    if ( !std::strcmp(ftype_shape, "rbeta" )) {
	       rshape = &rbeta;
	    } else if(!std::strcmp(ftype_shape,"const")) {
	       rshape = &rconst;
	    } else if(!std::strcmp(ftype_shape, "rlnorm")) {
	   	   rshape = &rlnorm;
	    } else if(!std::strcmp(ftype_shape,"rgamma")) {
	       rshape = &rgamma;
	    } else if(!std::strcmp(ftype_shape,"runif")) {
	       rshape = &runif;
	    } else {
	    	error(_("Unknown shape distribution type."));
	    }

	    rdist2_t rsize;
	    if( !std::strcmp( ftype_size, "rbinorm")) {
	    	simBivariate(R_args,rshape,dtype,label,isPerfect);
	    } else {
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
	}

	PutRNGstate();
	UNPROTECT(2);
	return;
}


/**
 * @brief Simulation of Spheroid system
 *        For constant size there is no perfect simulation!
 *
 *
 * @param R_param
 * @param R_cond
 * @return
 */
SEXP EllipsoidSystem(SEXP R_param, SEXP R_cond) {
  /* init */
  STGM::CSpheroidSystem *sp = InitSpheroidSystem(R_cond);

  /* simulate  */
  GetRNGstate();
  sp->simSpheroidSystem(R_param,R_cond);
  if(PL>10) Rprintf("Simulated %d spheroids: %p \n",sp->refObjects().size(), sp);

  SEXP R_spheroids = R_NilValue;
  PROTECT(R_spheroids = convert_R_EllipsoidSystem(sp));
  classgets(R_spheroids, getListElement( R_cond, "stype" ));
  PutRNGstate();

  _free_spheroids(sp);
  UNPROTECT(1);
  return R_spheroids;
}

SEXP IntersectSpheroidSystem(SEXP R_var, SEXP R_n, SEXP R_dz, SEXP R_intern, SEXP R_env, SEXP R_pl)
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
    error(_("Sphere system is missing a simulation box."));

  double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
  double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
  double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));
  STGM::CBox3 box(boxX,boxY,boxZ);

  SEXP R_mu = R_NilValue;
  PROTECT(R_mu = getAttrib(R_S,install("R_mu"))); ++nprotect;
  if(isNull(R_mu))
	error(_("Main orientation direction must not be 'Null'."));
  STGM::CVector3d mu(REAL(R_mu)[0],REAL(R_mu)[1],REAL(R_mu)[2]);

  const char* stype_str = GET_OBJECT_CLASS(R_S);
  STGM::CSpheroid::spheroid_type stype = STGM::CSpheroid::PROLATE;
   if( !std::strcmp("oblate",stype_str))
     stype = STGM::CSpheroid::OBLATE;

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
  STGM::CSpheroidSystem *sp = allocSpheroidSystem(box,lam,mu,stype,perfect);

  // do conversions
  sp->refObjects() = convert_C_Spheroids(R_S);

  if(PL>10) {
    Rprintf("Intersect with plane: %d , %p \n", sp->refObjects().size(), sp);
  }

  // intersect
  STGM::Intersectors<STGM::CSpheroid>::Type objects;
  sp->IntersectWithPlane(objects,plane,intern);

  SEXP R_ellipses;
  if(INTEGER(R_pl)[0] == 10) {									/* return short list of ellipses */
    PROTECT(R_ellipses = convert_R_Ellipses_trunc(objects)); ++nprotect;
  } else {
    PROTECT(R_ellipses = convert_R_Ellipses_all(objects)); ++nprotect;
  }

  SET_CLASS_NAME(R_ellipses,stype_str);

  _free_spheroids(sp);
  UNPROTECT(nprotect);
  return R_ellipses;
}


SEXP SimulateSpheroidsAndIntersect(SEXP R_param, SEXP R_cond)
{
  int nprotect = 0;
  /* init */
  STGM::CSpheroidSystem *sp = InitSpheroidSystem(R_cond);

  /* simulate  */
  sp->simSpheroidSystem(R_param,R_cond);
  if(PL>10) Rprintf("Simulated %d spheroids: %p \n",sp->refObjects().size(), sp);

  int intern = 0;
  SEXP R_intern = R_NilValue;
  PROTECT(R_intern = getListElement(R_cond,"intern")); ++nprotect;
  if(!isNull(R_intern))
	intern = INTEGER(AS_INTEGER(R_intern))[0];
  else { warning(_("Undefined argument `intern`. Set to zero.")); }

  double dz = 0.0;
  SEXP R_dz = R_NilValue;
  PROTECT(R_dz = AS_NUMERIC(getListElement(R_cond,"dz"))); ++nprotect;
  if(!isNull(R_dz))
	dz = REAL(R_dz)[0];
  else { warning(_("Undefined argument `dz`. Set to zero.")); }

  SEXP R_n = R_NilValue;
  PROTECT(R_n = AS_NUMERIC(getListElement(R_cond,"nsect"))); ++nprotect;
  if(isNull(R_n))
   error(_("Normal vector defining intersecting plane is `Null`."));
  STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
  STGM::CPlane plane(n,dz);

  STGM::Intersectors<STGM::CSpheroid>::Type objects;
  sp->IntersectWithPlane(objects,plane,intern);

  if(PL>10){
	  Rprintf("Plane normal to: [%f %f %f] \n", n[0],n[1],n[2]);
	  Rprintf("Number of intersections: %d at %f (intern=%d) \n",objects.size(),dz,intern);
  }

  SEXP R_ellipses;
  if(PL==10) {
	PROTECT(R_ellipses = convert_R_Ellipses_trunc(objects)); ++nprotect;
  } else {
	PROTECT(R_ellipses = convert_R_Ellipses_all(objects)); ++nprotect;
  }
  const char *stype_str=(sp->m_stype==0 ? "prolate" : "oblate");
  SET_CLASS_NAME(R_ellipses,stype_str);

  _free_spheroids(sp);
  UNPROTECT(nprotect);
  return R_ellipses;
}

SEXP convert_C2R_ellipses(STGM::Ellipses &ellipses) {
  int nProtected=0, dim=2, ncomps=8;
  int n = ellipses.size();

  SEXP names, R_resultlist;
  PROTECT(names = allocVector(STRSXP, ncomps));   ++nProtected;
  PROTECT(R_resultlist = allocVector(VECSXP,n));  ++nProtected;

  SET_STRING_ELT(names, 0, mkChar("id"));
  SET_STRING_ELT(names, 1, mkChar("center"));
  SET_STRING_ELT(names, 2, mkChar("ab"));
  SET_STRING_ELT(names, 3, mkChar("minor"));
  SET_STRING_ELT(names, 4, mkChar("major"));
  SET_STRING_ELT(names, 5, mkChar("A"));
  SET_STRING_ELT(names, 6, mkChar("phi"));
  SET_STRING_ELT(names, 7, mkChar("rot"));

  SEXP R_tmp,R_minor,R_major,R_A,R_center,R_ab;

  for(int i = 0; i < n; i++) {
      STGM::CEllipse2 & ellipse = ellipses[i];
      PROTECT(R_tmp = allocVector(VECSXP,ncomps));
      PROTECT(R_center = allocVector(REALSXP, dim));
      PROTECT(R_ab = allocVector(REALSXP, dim));
      PROTECT(R_minor = allocVector(REALSXP, dim));
      PROTECT(R_major = allocVector(REALSXP, dim));
      PROTECT(R_A = allocMatrix(REALSXP, dim,dim));

      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];

      REAL(R_minor)[0] = ellipse.minorAxis()[0];
      REAL(R_minor)[1] = ellipse.minorAxis()[1];

      REAL(R_major)[0] = ellipse.majorAxis()[0];
      REAL(R_major)[1] = ellipse.majorAxis()[1];

      REAL(R_ab)[0] = ellipse.a();
      REAL(R_ab)[1] = ellipse.b();

      COPY_C2R_MATRIX(ellipse.MatrixA(),R_A,dim);

      setAttrib(R_tmp, R_NamesSymbol, names);
      SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_ab);
      SET_VECTOR_ELT(R_tmp,3,R_minor);
      SET_VECTOR_ELT(R_tmp,4,R_major);
      SET_VECTOR_ELT(R_tmp,5,R_A);
      SET_VECTOR_ELT(R_tmp,6,ScalarReal(ellipse.phi()));
      SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.rot()));

      SET_VECTOR_ELT(R_resultlist,i,R_tmp);
      UNPROTECT(6);
  }

  UNPROTECT(nProtected);
  return R_resultlist;
}

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  SEXP R_tmp, R_ctr, R_u, R_ab, R_angles;

  SEXP R_label = R_NilValue;
  PROTECT(R_label = getAttrib(R_spheroids, install("label")));
  const char *label = (isNull(R_label) == TRUE ? "N" : translateChar(asChar(R_label)));

  int interior = 1;
  STGM::Spheroids spheroids;
  for(int i=0; i<length(R_spheroids); i++)
  {
      PROTECT(R_tmp = VECTOR_ELT(R_spheroids,i));

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

      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_tmp, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_tmp, "u")));
      PROTECT( R_ab     = AS_NUMERIC( getListElement( R_tmp, "acb")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_tmp, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
      STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

      spheroids.push_back(STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],REAL(R_ab)[2],u,
    		  REAL(R_angles)[0],REAL(R_angles)[1],INTEGER(AS_INTEGER( getListElement( R_tmp, "id")))[0],
			   label,interior));

      UNPROTECT(5);
  }

  UNPROTECT(1);
  return spheroids;
}

SEXP convert_R_Ellipses_trunc(STGM::Intersectors<STGM::CSpheroid>::Type &objects) {
  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp=R_NilValue;
  const char *nms[] = {"A", "C", "S", "phi", ""};

  double phi = 0.0;
  for(size_t k=0; k<objects.size(); ++k)   {
      STGM::CEllipse2 &ellipse = objects[k].getEllipse();
      PROTECT(R_tmp = mkNamed(VECSXP, nms));
      SET_VECTOR_ELT(R_tmp,0,ScalarReal(ellipse.a()));                  // major semi-axis (for both prolate/oblate)
      SET_VECTOR_ELT(R_tmp,1,ScalarReal(ellipse.b()));                  // minor semi-axis (for both prolate/oblate)
      SET_VECTOR_ELT(R_tmp,2,ScalarReal(ellipse.b()/ellipse.a()));      // shape

      phi = ellipse.phi();
      if(phi > M_PI_2) {
         if(phi <= M_PI) phi = M_PI-phi;
         else if(phi < 1.5*M_PI) phi = std::fmod(phi,M_PI);
         else phi = 2.0*M_PI-phi;
      }
      SET_VECTOR_ELT(R_tmp,3,ScalarReal(phi));
      SET_VECTOR_ELT(R_resultlist,k,R_tmp);
      UNPROTECT(1);
  }
  UNPROTECT(1);
  return(R_resultlist);
}

/**
 * @brief Convert ellipses to R objects
 * 		  Need angle between [0,2pi] for
 * 		  intersection plane relative to z axis
 * 		  but rgl uses relative to x axis!!
 * 		  For plotting: phi = phi + pi/2
 *
 * @param objects Intersection object
 * @return R ellipses
 */
SEXP convert_R_Ellipses_all(STGM::Intersectors<STGM::CSpheroid>::Type &objects) {
  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, objects.size()) );

  SEXP R_tmp, R_center, R_ab, R_A;
  const char *nms[] = {"id", "center", "A", "ab", "phi", "S", ""};

  for(size_t k=0; k<objects.size(); ++k)
  {
      STGM::CEllipse2 &ellipse = objects[k].getEllipse();

      PROTECT(R_tmp = mkNamed(VECSXP,nms));
      PROTECT(R_center = allocVector(REALSXP, 2));
      PROTECT(R_ab = allocVector(REALSXP, 2));
      PROTECT(R_A = allocMatrix(REALSXP,2,2));

      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];

      REAL(R_ab)[0] = ellipse.a();    // major semi-axis (for both prolate/oblate)
      REAL(R_ab)[1] = ellipse.b();	  // minor semi-axis (for both prolate/oblate)

      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
          REAL(R_A)[i + 2 *j] = ellipse.MatrixA()[i][j];

      SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_A);
      SET_VECTOR_ELT(R_tmp,3,R_ab);

      /* return angle between [0,2pi] */
      SET_VECTOR_ELT(R_tmp,4,ScalarReal(ellipse.phi()));
      SET_VECTOR_ELT(R_tmp,5,ScalarReal(ellipse.b()/ellipse.a()));
      SET_VECTOR_ELT(R_resultlist,k,R_tmp);
      UNPROTECT(4);
  }

  UNPROTECT(1);
  return(R_resultlist);
}

STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  SEXP R_ctr, R_u, R_acb, R_angles;
  PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_spheroid, "center")));
  PROTECT( R_u      = AS_NUMERIC( getListElement( R_spheroid, "u")));
  PROTECT( R_acb    = AS_NUMERIC( getListElement( R_spheroid, "acb")));
  PROTECT( R_angles = AS_NUMERIC( getListElement( R_spheroid, "angles")));

  STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
  STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

  int interior = 1;
  const char *label = "N";
  if(!isNull(getAttrib(R_spheroid, install("label"))))
    label = translateChar(asChar(getAttrib(R_spheroid, install("label"))));
  else {
	  label = "N";
	  warning(_("Undefined attribute `label`. Set to 'N'."));
  }
  if(!isNull(getAttrib(R_spheroid, install("interior"))))
    interior = asLogical(getAttrib(R_spheroid, install("interior")));
  else {
	interior = 0;
	warning(_("Undefined attribute `interior`. Set to zero."));
  }

  UNPROTECT(4);
  return STGM::CSpheroid(ctr,REAL(R_acb)[0],REAL(R_acb)[1],REAL(R_acb)[2],u,
		   	   	   	   	 REAL(R_angles)[0],REAL(R_angles)[1],
						 INTEGER(AS_INTEGER(getListElement( R_spheroid, "id")))[0],
						 label,interior);
}


/**
 * Needed for digitization only
 */
/*
STGM::Ellipses convert_C_Ellipses(SEXP R_ellipses)
{
  STGM::Ellipses ellipses;
  SEXP R_tmp, R_ctr,R_ab, R_minor, R_major;

  for(int i=0; i<length(R_ellipses); i++) {
     PROTECT(R_tmp   = VECTOR_ELT(R_ellipses,i));
     PROTECT(R_ctr   = AS_NUMERIC( getListElement( R_tmp, "center")));
     PROTECT(R_minor = AS_NUMERIC( getListElement( R_tmp, "minor")));
     PROTECT(R_major = AS_NUMERIC( getListElement( R_tmp, "major")));
     PROTECT(R_ab    = AS_NUMERIC( getListElement( R_tmp, "ab")));

    //BUG: Constructor with matrix A has a bug to determine the correct angle phi

     STGM::CPoint2d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1]);
     STGM::CPoint2d minorAxis(REAL(R_minor)[0],REAL(R_minor)[1]);
     STGM::CPoint2d majorAxis(REAL(R_major)[0],REAL(R_major)[1]);

     // m_A is re-computed
     ellipses.push_back(STGM::CEllipse2(ctr,majorAxis,minorAxis,
    		 	 	 	REAL(R_ab)[0],REAL(R_ab)[1],
						INTEGER(AS_INTEGER( getListElement( R_tmp, "id")))[0],
						REAL(AS_NUMERIC(getListElement( R_tmp, "A")))[0]));

     UNPROTECT(5);
  }

  return ellipses;
}
*/


SEXP convert_R_EllipsoidSystem( STGM::CSpheroidSystem *S) {
  int nProtected=0, ncomps=6, dim=3;

  STGM::CBox3 &box = S->box();
  STGM::Spheroids &spheroids = S->refObjects();

  SEXP R_resultlist;
  PROTECT(R_resultlist = allocVector(VECSXP, spheroids.size()) ); ++nProtected;

  SEXP names, R_center, R_u, R_acb, R_tmp, R_angles, R_rotM;
  // get lateral bounding planes
  const STGM::LateralPlanes &planes = box.getLateralPlanes();

  for(size_t k=0;k<spheroids.size();k++)
  {
    STGM::CSpheroid &spheroid = spheroids[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, dim));
    PROTECT(R_u = allocVector(REALSXP, dim));
    PROTECT(R_acb = allocVector(REALSXP,3));
    PROTECT(R_angles = allocVector(REALSXP,2));
    PROTECT(R_rotM = allocMatrix(REALSXP,3,3));

    // projection
    //STGM::CEllipse2 ellipse = spheroid.spheroidProjection();
    STGM::CEllipse2 ellipse = spheroid.delamProjection();

    // check intersection with box (surface crack or inner damage)
    STGM::Intersector<STGM::CSpheroid> intersector(spheroid , box.m_size );
    Rboolean interior = (Rboolean) TRUE;
    for(size_t j=0; j<planes.size() ; ++j) {
         if( intersector(planes[j])) {
           interior = (Rboolean) FALSE;
           break;
         }
    }

    const STGM::CVector3d &m_center = spheroid.center();
    REAL(R_center)[0]=m_center[0];
    REAL(R_center)[1]=m_center[1];
    REAL(R_center)[2]=m_center[2];

    const STGM::CVector3d &m_u = spheroid.u();
    REAL(R_u)[0]=m_u[0];
    REAL(R_u)[1]=m_u[1];
    REAL(R_u)[2]=m_u[2];

    REAL(R_acb)[0]=spheroid.a();
    REAL(R_acb)[1]=spheroid.c();
    REAL(R_acb)[2]=spheroid.b();    // major semi-axis

    REAL(R_angles)[0]=spheroid.theta();
    REAL(R_angles)[1]=spheroid.phi();

    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
        REAL(R_rotM)[i + dim *j] = (spheroid.rotationMatrix())[i][j];

    //STGM::CMatrix3d &M = spheroid.rotationMatrix();
    //COPY_C2R_MATRIX(M,R_rotM,dim);

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(spheroid.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,R_u);
    SET_VECTOR_ELT(R_tmp,3,R_acb);
    SET_VECTOR_ELT(R_tmp,4,R_angles);
    SET_VECTOR_ELT(R_tmp,5,R_rotM);

    PROTECT(names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(names, 0, mkChar("id"));
    SET_STRING_ELT(names, 1, mkChar("center"));
    SET_STRING_ELT(names, 2, mkChar("u"));
    SET_STRING_ELT(names, 3, mkChar("acb"));
    SET_STRING_ELT(names, 4, mkChar("angles"));
    SET_STRING_ELT(names, 5, mkChar("rotM"));
    setAttrib(R_tmp, R_NamesSymbol, names);
    setAttrib(R_tmp, install("label"), mkString(spheroid.label()) );
    setAttrib(R_tmp, install("interior"), ScalarLogical(interior));
    setAttrib(R_tmp, install("area"), ScalarReal(ellipse.area()));

    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(7);
  }

  UNPROTECT(nProtected);
  return(R_resultlist);
}

/**
 * @brief simulation spheroid system but not perfect
 *        using user defined random joint distribution
 *
 * @param d R calling data
 *
 * */
void STGM::CSpheroidSystem::simJoint(SEXP R_call, SEXP R_rho, const char *label) {
     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
       num = rpois(m_box.volume()*m_lam);
       ++nTry;
     }
     m_spheroids.reserve(num);

     double m1 = m_box.m_size[0] +(m_box.m_center[0]-m_box.m_extent[0]),
            m2 = m_box.m_size[1] +(m_box.m_center[1]-m_box.m_extent[1]),
            m3 = m_box.m_size[2] +(m_box.m_center[2]-m_box.m_extent[2]);

     double *v=0, a=0, b=0, c=0, theta=0, phi=0;

     CVector3d u;
     SEXP Reval = R_NilValue;

     int err = 0;
     for (size_t niter=0; niter<num; niter++)
     {
         Reval = R_tryEval(R_call,R_rho,&err);
         if(!err) {
            a=REAL(getListElement(Reval,"a"))[0];
            b=REAL(getListElement(Reval,"b"))[0];
            c=REAL(getListElement(Reval,"c"))[0];
            theta=REAL(getListElement(Reval,"theta"))[0];
            phi=REAL(getListElement(Reval,"phi"))[0];

            v=REAL(getListElement(Reval,"u"));
            u[0]=v[0]; u[1]=v[1]; u[2]=v[2];

            if(m_stype==CSpheroid::OBLATE)
              std::swap(a,b);

            STGM::CVector3d center(runif(0.0,1.0)*m1,runif(0.0,1.0)*m2, runif(0.0,1.0)*m3);
            m_spheroids.push_back( STGM::CSpheroid(center,a,c,b,u,theta,phi,m_spheroids.size()+1,label) );
         } else
            error(_("simJoint(): R `try` error in user defined distribution function."));
     }

}

/**
 * @brief Bivariate ellipsoid-size-shape distribution,
 *        the major semiaxis is logN distributed,
 *        and possibly shorter semiaxes of unequal lengths
 *
 *
 * @param d data of R call (no R call used here)
 */
void STGM::CSpheroidSystem::simBivariate(SEXP R_args, rdist2_t rshape, STGM::CSpheroid::direction_type dtype,
		                                    const char *label, int perfect)
{
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

      // shape distribution parameters for both shorter semi-axes
      // defaults are s1=s2=1.0 in R
      double s1 = REAL(VECTOR_ELT(VECTOR_ELT(R_args,1),0))[0];
      double s2 = LENGTH(VECTOR_ELT( R_args, 1)) > 1 ? REAL(VECTOR_ELT(VECTOR_ELT( R_args, 1),1))[0] : 0;

      if(PL>10) {
         Rprintf("Spheroids (exact) simulation: \n");
         Rprintf("\t size distribution: `rbinorm`  with mx=%f, sdx=%f, my=%f, sdy=%f, rho=%f, kappa=%f \n",mx,sdx,my,sdy,rho,kappa);
         Rprintf("\t directional distribution: %d  \n", dtype);
         Rprintf("\t cum sum of probabilities: %f, %f, %f, %f \n",p[0],p[1],p[2],p[3]);
         Rprintf("\t set label: %s to character: \n",label);
         Rprintf("\t Shape parameters for shorter semi-axes: s1 = %f, s2 = %f \n", s1, s2);
         Rprintf("\n\n");
      }

      int nTry=0;
      while(num==0 && nTry<MAX_ITER) {
         num = rpois(mu*m_lam);
         ++nTry;
      }
      m_spheroids.reserve(num);

      CVector3d u;
      double x=0,y=0,r=0,
    		 a=0,c=0,b=0,    						/* two shorter semiaxes a,c and major semiaxis b */
    		 s=1.0,phi=0,theta=0;

      if(perfect) {

    	  for (size_t niter=0; niter<num; niter++)
    	  {
    		  rbinorm_exact(p,mx,sdx,my,sdy,rho,x,y);
    		  s=1.0/(1.0+std::exp(-y));
    		  b=r=std::exp(x); 						/* b = r for exact simulation*/
    		  a=b*s;

    		  c=a*rshape(s1,s2);					/* either constant factor or random 2nd. shorter axis */
			  switch(dtype) {					    /* sample orientation */
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

			  STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
			                         runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
			                         runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));

			  /* b = r */
			  m_spheroids.push_back( STGM::CSpheroid(center,a,c,b,u,theta,phi,niter+1,label) );
    	  }

      } else {

    	  for (size_t niter=0; niter<num; niter++)
    	  {
    		  rbinorm(mx,sdx,my,sdy,rho,x,y);
    		  s=1.0/(1.0+std::exp(-y));
    		  b=std::exp(x); 						/* b = r for exact simulation*/
    		  a=b*s;
    		  if(m_stype==CSpheroid::OBLATE)
				  std::swap(a,b);
			  /* either constant or random to allow for general ellipsoids */
			  c=a*rshape(s1,s2);

			  /* sample orientation */
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

			  STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
								   	 runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
									 runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

			  m_spheroids.push_back( STGM::CSpheroid(center,a,c,b,u,theta,phi,niter+1,label) );

    	  }
      }

}



/**
 * @brief Perfect simulation for (non constant) size distributions
 *        with independent shape and orientation distributions
 *
 * @param d    R call data
 */
void STGM::CSpheroidSystem::simUnivar(SEXP R_args, rdist2_t rsize, rdist2_t rshape,
							STGM::CSpheroid::direction_type dtype, const char *label)
{
     int nTry=0;
     while(num==0 && nTry<MAX_ITER) {
        num = rpois(m_box.volume()*m_lam);
        ++nTry;
     }
     m_spheroids.reserve(num);

     /// size
	 double p1=REAL(VECTOR_ELT(VECTOR_ELT( R_args, 0),0))[0];
	 double p2=LENGTH(VECTOR_ELT( R_args, 0)) > 1 ? REAL(VECTOR_ELT(VECTOR_ELT( R_args, 0),1))[0] : 0;

	 // shape
	 double s1=REAL(VECTOR_ELT(VECTOR_ELT(R_args,1),0))[0];
	 double s2=LENGTH(VECTOR_ELT( R_args, 1)) > 1 ? REAL(VECTOR_ELT(VECTOR_ELT( R_args, 1),1))[0] : 0;

     // direction
     double kappa = REAL(VECTOR_ELT(VECTOR_ELT(R_args,2),0))[0];

     if(PL>10) {
         Rprintf("Run spheroids  simulation... \n");
         Rprintf("\t size distribution: %f %f \n", p1,p2);
         Rprintf("\t shape parameters: %f %f\n", s1, s2);
         Rprintf("\t directional distribution: %d  with %f \n", dtype, kappa);
         Rprintf("\t set label: %s to character: \n",label);
     }

     /* loop over all */
     CVector3d u;
     double a=0, b=0, theta=0, phi=0;
     for (size_t niter=0; niter<num; niter++) {
    	 b = rsize(p1,p2);      			/* major semi-axis */
         a = b * rshape(s1,s2); 			/* minor semi-axis, shape factor, constant or random */
         if(m_stype==CSpheroid::OBLATE)
           std::swap(a,b);					/* just swap a with b */

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

        STGM::CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
                               runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
                               runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

        m_spheroids.push_back( STGM::CSpheroid(center,a,a,b,u,theta,phi,niter+1,label) );
     }
}

/**
 * @brief Intersect only objects whose center is inside the window
 *
 * @param objects
 * @param plane
 * @param intern
 */
void STGM::CSpheroidSystem::IntersectWithPlane(STGM::Intersectors<STGM::CSpheroid>::Type &objects,
											   STGM::CPlane &plane, int intern)
{
  int i=0,j=0;
  switch(plane.idx()) {
      case 0: i=1; j=2; break; // YZ
      case 1: i=0; j=2; break; // XZ
      case 2: i=0; j=1; break; // XY
  }

  // assume left-down corner is origin of box
  if(intern) {
      CWindow win(m_box.m_size[i],m_box.m_size[j]);
      for(size_t i=0; i<m_spheroids.size(); ++i) {
           STGM::Intersector<STGM::CSpheroid> intersector( m_spheroids[i], plane, m_box.m_size);
           if(intersector.FindIntersection()) {
               if(intersector.getEllipse().isInWindow(win))
                 objects.push_back( intersector );
           }
       }
  } else {
      for(size_t i=0; i<m_spheroids.size(); ++i) {
          STGM::Intersector<STGM::CSpheroid> intersector( m_spheroids[i], plane, m_box.m_size);
          if(intersector.FindIntersection()) {
            objects.push_back( intersector );
          }
      }
  }
}

/**
 * Test whether simulation box is intersected
 */
SEXP UpdateIntersections(SEXP R_S, SEXP R_box) {
    int nProtected=0;

    SEXP R_BoxX, R_BoxY, R_BoxZ;
    PROTECT( R_BoxX = AS_NUMERIC( getListElement( R_box, "xrange" ) ) ); ++nProtected;
    PROTECT( R_BoxY = AS_NUMERIC( getListElement( R_box, "yrange" ) ) ); ++nProtected;
    PROTECT( R_BoxZ = AS_NUMERIC( getListElement( R_box, "zrange" ) ) ); ++nProtected;

    STGM::CBox3 box(REAL(R_BoxX)[1],REAL(R_BoxY)[1],REAL(R_BoxZ)[1]);
    const std::vector<STGM::CPlane> &planes = box.getPlanes();

    SEXP R_ret = R_NilValue;
    PROTECT(R_ret = allocVector(INTSXP,length(R_S)));  ++nProtected;

    const char *name = GET_OBJECT_CLASS(R_S);
    if( !std::strcmp(name, "prolate" ) || !std::strcmp(name, "oblate" ) || !std::strcmp(name, "spheroid" )) {
        for(int k=0;k<length(R_S);k++) {
                STGM::CSpheroid sp = convert_C_Spheroid(VECTOR_ELT(R_S,k));
                STGM::Intersector<STGM::CSpheroid> intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!std::strcmp(name, "cylinder" )) {
        for(int k=0;k<length(R_S);k++) {
                STGM::CCylinder sp = convert_C_Cylinder(VECTOR_ELT(R_S,k));
                STGM::Intersector<STGM::CCylinder> intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!std::strcmp(name, "sphere" )) {
      for(int k=0;k<length(R_S);k++) {
                STGM::CSphere sp = convert_C_Sphere(VECTOR_ELT(R_S,k));
                STGM::Intersector<STGM::CSphere> intersector(sp , box.m_size );
                INTEGER(R_ret)[k] = intersector.TestBoxIntersection(planes);
      }
    } else {
        error(_("Unknown class object."));
    }

    UNPROTECT(nProtected);
    return R_ret;
}


SEXP DigitizeProfiles(SEXP R_S, SEXP R_cond, SEXP R_delta)
{
	int nProtected = 0;
	SEXP R_ret = R_NilValue;
	PROTECT(R_ret = allocVector(INTSXP,length(R_S)));  ++nProtected;
	const char *name = GET_OBJECT_CLASS(R_S);

	if( !std::strcmp(name, "prolate" ) || !std::strcmp(name, "oblate" ) || !std::strcmp(name, "spheroid" )) {
		for(int k=0;k<length(R_S);k++) {

		}
	} else if(!std::strcmp(name, "cylinder" )) {
		for(int k=0;k<length(R_S);k++) {

		}
	} else if(!std::strcmp(name, "sphere" )) {
	  for(int k=0;k<length(R_S);k++) {

	  }
	} else {
		error(_("Unknown class object."));
	}

	 UNPROTECT(nProtected);
	 return R_ret;

  /*
   *
   checkPtr(ext, spheroid_type_tag);
   STGM::CSpheroidSystem *sp = static_cast<STGM::CSpheroidSystem *>(getExternalPtr(ext));
   STGM::CVector3d n(REAL(R_n)[0],REAL(R_n)[1],REAL(R_n)[2]);
   STGM::CPlane plane( n , asReal(R_z));

   if(PL>10) Rprintf("Intersect with plane: %d \n", sp->refObjects().size());
   STGM::Intersectors<STGM::CSpheroid>::Type objects;
   sp->IntersectWithPlane(objects,plane,0);

   if(PL>10) Rprintf("done: %d \n", objects.size());
   int nPix = (int) sp->box().m_size[0]/REAL(R_delta)[0];

   SEXP R_W = R_NilValue;
   PROTECT(R_W = allocMatrix(INTSXP,nPix,nPix));

   if(PL>10) Rprintf("Digitize: nPix: %d, delta: %f \n",nPix,REAL(R_delta)[0]);
   STGM::digitize<STGM::CSpheroid>(objects,INTEGER(R_W),nPix,asReal(R_delta));

   UNPROTECT(1);
   return R_W;
*/

 }
