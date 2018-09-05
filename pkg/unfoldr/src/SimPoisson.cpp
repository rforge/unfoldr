/**
 *  @file SimEllipsoid.cpp
 *
 *  @date: 10.01.2014
 *  @author: M. Baaske
 */

#define MAX_ITER 100

#include "directions.h"
#include "SimPoisson.h"

//static locals
static int PL = 0;

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

#define SET_REAL_VECTOR(Rv, vec)  			\
do {										\
  int _i;		                            \
  for(_i = 0; _i < LENGTH( (Rv) ); _i++) {	\
    REAL((Rv))[_i] = (vec)[_i];				\
  }											\
} while(0)


/* auxiliary functions */
STGM::CBox3 setBox(SEXP R_box);
void setWindow(SEXP R_win, STGM::CBox3 &box, STGM::CPlane &plane);


using namespace std;

namespace STGM {

 template<typename T>
 void IntersectWithPlane(CPoissonSystem<T> &sp, typename Intersectors<T>::Type &intersected, SEXP R_cond);

}



/**
 * @brief Simulation of Poisson system
 *        For constant size there is no perfect simulation!
 */
SEXP PoissonSystem(SEXP R_param, SEXP R_cond) {
	int nprotect = 0;
	SEXP R_ret = R_NilValue;
	const char* type_str = CHAR( STRING_ELT( getListElement( R_cond, "type" ), 0 ));
	const char* profiles = CHAR( STRING_ELT( getListElement( R_cond, "intersect"), 0 ));

	SEXP R_box = R_NilValue;
	PROTECT( R_box  = getListElement( R_cond, "box"));	++nprotect;
	STGM::CBox3 box = setBox(R_box);

	// set print level
	PL = INTEGER(AS_INTEGER(getListElement( R_cond,"pl")))[0];

	SEXP R_exact;
	PROTECT(R_exact = getListElement( R_cond,"perfect")); ++nprotect;
	int perfect = INTEGER(AS_INTEGER(R_exact))[0];

	SEXP R_lam;
	PROTECT(R_lam = getListElement( R_cond, "lam")); ++nprotect;
	double lam = REAL(AS_NUMERIC(R_lam))[0];

	/* set up spheroid system */
	SEXP R_mu;
	PROTECT(R_mu = getListElement( R_cond, "mu")); ++nprotect;
	STGM::CVector3d maxis(REAL(AS_NUMERIC(R_mu)));

	if( !std::strcmp(type_str, "prolate" ) ||
		!std::strcmp(type_str, "oblate" ))
	{
		STGM::CPoissonSystem<STGM::CSpheroid> sp(box,lam,maxis,type_str,perfect);
		sp.simSystem(R_param,R_cond);

		if(!std::strcmp(profiles, "full") ||
		   !std::strcmp(profiles, "only") ||
		   !std::strcmp(profiles, "digi"))
		{
			STGM::Intersectors<STGM::CSpheroid>::Type intersected;
			STGM::IntersectWithPlane<STGM::CSpheroid>(sp, intersected, R_cond);

			if(!std::strcmp(profiles, "full"))
			{
				const char *nms[] = {"S", "sp", "W", ""};
				PROTECT(R_ret = mkNamed(VECSXP, nms));
				SET_VECTOR_ELT(R_ret,0,convert_R_Ellipsoids(sp));
				setAttrib( VECTOR_ELT(R_ret,0), install("mu"), R_mu);
				setAttrib( VECTOR_ELT(R_ret,0), install("lam"), R_lam);
				setAttrib( VECTOR_ELT(R_ret,0), install("box"), R_box);
				setAttrib( VECTOR_ELT(R_ret,0), install("perfect"), R_exact);

				/* intersection */
				SET_VECTOR_ELT(R_ret,1,convert_R_Ellipses(intersected, box));

				/* digitization */
				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
										  (int) std::floor(win.m_size[1]/delta));

				SEXP R_tmp = R_NilValue;
				PROTECT(R_tmp = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_tmp),win.m_low,nPix,delta);
				digitizer.start<STGM::CSpheroid>(intersected);

				SET_VECTOR_ELT(R_ret,2,R_tmp);
				UNPROTECT(1);

			}  else if(!std::strcmp(profiles, "digi" )) {			/* only digitized image */

				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
							  		      (int) std::floor(win.m_size[1]/delta));

				PROTECT(R_ret = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_ret),win.m_low,nPix,delta);
				digitizer.start<STGM::CSpheroid>(intersected);

			} else {
				PROTECT(R_ret = convert_R_Ellipses(intersected, box));
			}

		} else if(!std::strcmp(profiles, "original" )) {

			PROTECT(R_ret = convert_R_Ellipsoids(sp) );
			setAttrib( R_ret, install("mu"), R_mu);
			setAttrib( R_ret, install("lam"), R_lam);
			setAttrib( R_ret, install("box"), R_box);
			setAttrib( R_ret, install("perfect"), R_exact);

		}

	} else if(!std::strcmp(type_str, "cylinders" )) {

		/* init */
		STGM::CPoissonSystem<STGM::CCylinder> sp(box,lam,maxis,type_str,perfect);
		sp.simSystem(R_param,R_cond);

		if(!std::strcmp(profiles, "full") ||
		   !std::strcmp(profiles, "only") ||
		   !std::strcmp(profiles, "digi"))
		{
			STGM::Intersectors<STGM::CCylinder>::Type intersected;
			STGM::IntersectWithPlane<STGM::CCylinder>(sp,intersected,R_cond);

			if(!std::strcmp(profiles, "full"))
			{
				const char *nms[] = {"S", "sp", "W", ""};
				PROTECT(R_ret = mkNamed(VECSXP, nms));
				SET_VECTOR_ELT(R_ret,0,convert_R_Cylinders(sp));
				setAttrib( VECTOR_ELT(R_ret,0), install("mu"), R_mu);
				setAttrib( VECTOR_ELT(R_ret,0), install("lam"), R_lam);
				setAttrib( VECTOR_ELT(R_ret,0), install("box"), R_box);
				setAttrib( VECTOR_ELT(R_ret,0), install("perfect"), R_exact);

				SET_VECTOR_ELT(R_ret,1,convert_R_CylinderIntersections(intersected, box));

				/* digitization */
				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
										  (int) std::floor(win.m_size[1]/delta));

				SEXP R_tmp = R_NilValue;
				PROTECT(R_tmp = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_tmp),win.m_low,nPix,delta);
				digitizer.start<STGM::CCylinder>(intersected);

				SET_VECTOR_ELT(R_ret,2,R_tmp);
				UNPROTECT(1);

			}  else if(!std::strcmp(profiles, "digi" )) {			/* only digitized image */

				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
										  (int) std::floor(win.m_size[1]/delta));

				PROTECT(R_ret = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_ret),win.m_low,nPix,delta);
				digitizer.start<STGM::CCylinder>(intersected);


			} else {
				PROTECT(R_ret = convert_R_CylinderIntersections(intersected, box));
			}

		}  else if(!std::strcmp(profiles, "original" )) {
			PROTECT(R_ret = convert_R_Cylinders(sp));
			setAttrib( R_ret, install("mu"), R_mu);
			setAttrib( R_ret, install("lam"), R_lam);
			setAttrib( R_ret, install("box"), R_box);
			setAttrib( R_ret, install("perfect"), R_exact);

		}

	} else if(!std::strcmp(type_str, "spheres" )) {

		STGM::CPoissonSystem<STGM::CSphere> sp(box,lam,maxis,type_str,perfect);
		sp.simSystem( R_param, R_cond);

		if(!std::strcmp(profiles, "full") ||
		   !std::strcmp(profiles, "only") ||
		   !std::strcmp(profiles, "digi"))
		{
			STGM::Intersectors<STGM::CSphere>::Type intersected;
			STGM::IntersectWithPlane<STGM::CSphere>(sp,intersected,R_cond);

			if(!std::strcmp(profiles, "full"))
			{
				const char *nms[] = {"S", "sp", "W", ""};
				PROTECT(R_ret = mkNamed(VECSXP, nms));
				SET_VECTOR_ELT(R_ret,0,convert_R_Spheres(sp));
				setAttrib( VECTOR_ELT(R_ret,0), install("mu"), R_mu);
				setAttrib( VECTOR_ELT(R_ret,0), install("lam"), R_lam);
				setAttrib( VECTOR_ELT(R_ret,0), install("box"), R_box);
				setAttrib( VECTOR_ELT(R_ret,0), install("perfect"), R_exact);

				SET_VECTOR_ELT(R_ret,1,convert_R_Circles(intersected, box));


				/* digitization */
				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
										  (int) std::floor(win.m_size[1]/delta));

				SEXP R_tmp = R_NilValue;
				PROTECT(R_tmp = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_tmp),win.m_low,nPix,delta);
				digitizer.start<STGM::CSphere>(intersected);

				SET_VECTOR_ELT(R_ret,2,R_tmp);
				UNPROTECT(1);


			} else if(!std::strcmp(profiles, "digi" )) {			/* only digitized image */

				STGM::CPlane & plane = intersected[0].getPlane();
				int i=0, j=0;
				plane.getPlaneIdx(i,j);
				STGM::CWindow win(box.m_size[i],box.m_size[j]);

				double delta = REAL(AS_NUMERIC(getListElement( R_cond, "delta")))[0];
				STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
										  (int) std::floor(win.m_size[1]/delta));

				PROTECT(R_ret = allocMatrix(INTSXP,nPix[0],nPix[1]));
				STGM::CDigitizer digitizer(INTEGER(R_ret),win.m_low,nPix,delta);
				digitizer.start<STGM::CSphere>(intersected);


			} else {
				/* return full circle object */
				PROTECT(R_ret = convert_R_Circles(intersected, box));
			}

		} else {
			PROTECT(R_ret = convert_R_Spheres(sp));
			setAttrib( R_ret, install("mu"), R_mu);
			setAttrib( R_ret, install("lam"), R_lam);
			setAttrib( R_ret, install("box"), R_box);
			setAttrib( R_ret, install("perfect"), R_exact);
	   }

	} else {
	   error(_("Unknown object type."));
	}

	UNPROTECT(nprotect+1);
    return R_ret;
}

SEXP IntersectPoissonSystem(SEXP R_var, SEXP R_cond, SEXP R_env)
{
  int nprotect = 0;
  SEXP R_S = R_NilValue;
  PROTECT(R_S = getVar(R_var,R_env)); ++nprotect;

  /* read all from attributes */
  SEXP R_box = R_NilValue;
  PROTECT(R_box = getAttrib(R_S, install("box"))); ++nprotect;
  STGM::CBox3 box = setBox(R_box);

  SEXP R_mu = R_NilValue;
  PROTECT(R_mu = getAttrib(R_S,install("mu")));  ++nprotect;
  if(isNull(R_mu))
	error(_("`mu` main direction must be provided as an attribute."));
  STGM::CVector3d maxis(REAL(AS_NUMERIC(R_mu)));

  SEXP R_lam = R_NilValue;
  PROTECT(R_lam = getAttrib(R_S,install("lam"))); ++nprotect;
  if(isNull(R_lam))
    error(_("`lam` must be provided as an attribute."));
  double lam = REAL(AS_NUMERIC(R_lam))[0];

  SEXP R_perfect = R_NilValue;
  PROTECT(R_perfect = getAttrib(R_S,install("perfect"))); ++nprotect;
  if(isNull(R_perfect))
   error(_("`perfect` must be provided as an attribute."));
  int perfect = INTEGER(AS_INTEGER(R_perfect))[0];

  /* get class */
  SEXP Rclass = R_NilValue;
  PROTECT(Rclass = getAttrib( R_S, R_ClassSymbol));
  const char* type_str = CHAR( STRING_ELT(Rclass,0));
  UNPROTECT(1);

  /* global print level */
  PL = INTEGER(AS_INTEGER(getListElement( R_cond,"pl")))[0];

  /* return value */
  SEXP R_ret = R_NilValue;
  if( !std::strcmp(type_str, "prolate" ) ||
  		!std::strcmp(type_str, "oblate" ))
  	{
  		STGM::CPoissonSystem<STGM::CSpheroid> sp(box,lam,maxis,type_str,perfect);
     	// do conversions
	    sp.refObjects() = convert_C_Spheroids(R_S);

	    // intersect
	    STGM::Intersectors<STGM::CSpheroid>::Type intersected;
	    STGM::IntersectWithPlane<STGM::CSpheroid>(sp, intersected, R_cond);

	   	PROTECT(R_ret = convert_R_Ellipses(intersected, box));

  } else if(!std::strcmp(type_str, "cylinders" )) {

		STGM::CPoissonSystem<STGM::CCylinder> sp(box,lam,maxis,type_str,perfect);
		sp.refObjects() = convert_C_Cylinders(R_S);

		// intersect
		STGM::Intersectors<STGM::CCylinder>::Type intersected;
		STGM::IntersectWithPlane<STGM::CCylinder>(sp, intersected, R_cond);

		PROTECT(R_ret = convert_R_CylinderIntersections(intersected, box));

  } else if(!std::strcmp(type_str, "spheres" )) {

	  STGM::CPoissonSystem<STGM::CSphere> sp(box,lam,maxis,type_str,perfect);
	  sp.refObjects() = convert_C_Spheres(R_S);
	  // intersect
	  STGM::Intersectors<STGM::CSphere>::Type intersected;
	  STGM::IntersectWithPlane<STGM::CSphere>(sp, intersected, R_cond);

	  PROTECT(R_ret = convert_R_Circles(intersected, box));

  } else {
	   error(_("Unknown simulation object type."));
  }

  UNPROTECT(nprotect+1);
  return R_ret;

}

STGM::CBox3 setBox(SEXP R_box) {
	if(isNull(R_box)){
		error(_("Undefined simulation box."));
	}
	double *boxX = NUMERIC_POINTER( getListElement( R_box, "xrange"));
	double *boxY = NUMERIC_POINTER( getListElement( R_box, "yrange"));
	double *boxZ = NUMERIC_POINTER( getListElement( R_box, "zrange"));

	return STGM::CBox3(boxX,boxY,boxZ);
}

/**
 * Test whether simulation box is intersected
 */
SEXP UpdateIntersections(SEXP R_var, SEXP R_env)
{
    SEXP R_S = R_NilValue;
    PROTECT(R_S = getVar(R_var,R_env));

    SEXP R_box = R_NilValue;
    PROTECT(R_box = getAttrib(R_S, install("box")));
    STGM::CBox3 box = setBox(R_box);
    UNPROTECT(1);

    const std::vector<STGM::CPlane> &planes = box.getPlanes();

    SEXP R_ret = R_NilValue;
    PROTECT(R_ret = allocVector(INTSXP,length(R_S)));
    int *ret = INTEGER(R_ret);

    /* get class */
    SEXP Rclass = R_NilValue;
    PROTECT(Rclass = getAttrib( R_S, R_ClassSymbol));
    const char* name = CHAR( STRING_ELT(Rclass,0));
    UNPROTECT(1);

    if( !std::strcmp(name, "prolate")
     || !std::strcmp(name, "oblate" )) {

    	for(int k=0;k<length(R_S);k++) {
        	STGM::CSpheroid sp = convert_C_Spheroid(VECTOR_ELT(R_S,k));
        	STGM::Intersector<STGM::CSpheroid> intersector(sp, box.m_size );
            ret[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!std::strcmp(name, "cylinders" )) {

    	for(int k=0;k<length(R_S);k++) {
        	STGM::CCylinder sp = convert_C_Cylinder(VECTOR_ELT(R_S,k));
        	STGM::Intersector<STGM::CCylinder> intersector(sp, box.m_size );
            ret[k] = intersector.TestBoxIntersection(planes);
        }
    } else if(!std::strcmp(name, "spheres" )) {

      for(int k=0;k<length(R_S);k++) {
    	  STGM::CSphere sp = convert_C_Sphere(VECTOR_ELT(R_S,k));
    	  STGM::Intersector<STGM::CSphere> intersector(sp, box.m_size );
          ret[k] = intersector.TestBoxIntersection(planes);
      }
    } else {
        error(_("Unknown class object."));
    }

    UNPROTECT(2);
    return R_ret;
}


/**
 * Digitize intersections
 *  @param
 */
SEXP DigitizeProfiles(SEXP R_var, SEXP R_delta, SEXP R_win, SEXP R_env)
{
	int nprotect = 0;
	SEXP R_S = R_NilValue;
	PROTECT(R_S = getVar(R_var,R_env)); ++nprotect;			// section profiles

	/*  get the window if not provided:
	 *  this has the correct dimension of the original box
	 *  and corresponds to the intersecting plane  (normal vector) */
	if(isNull(R_win)){
		PROTECT(R_win = getAttrib(R_S, install("win")));
		++nprotect;
		if(isNull(R_win))
		 error(_("Intersection window must be given as an attribute."));
	}

	/** intersection window */
	STGM::CWindow win(REAL(VECTOR_ELT(R_win,0)),REAL(VECTOR_ELT(R_win,1)));

	/** spacing */
	double delta = REAL(AS_NUMERIC(R_delta))[0];
	STGM::CVector<int,2> nPix((int) std::floor(win.m_size[0]/delta),
							  (int) std::floor(win.m_size[1]/delta));

	if(PL>10) {
	 Rprintf("Digitize intersection profiles with resolution (%d x %d) and spacing delta by %f \n",nPix[0],nPix[1], delta);
	}

	/* alloc binary return  matrix */
	SEXP R_w = R_NilValue;
	PROTECT(R_w = allocMatrix(INTSXP,nPix[0],nPix[1])); ++nprotect;

	/* init digitizer, pass lower left corner in order to move to objects relative to [0,0] */
	STGM::CDigitizer digitizer(INTEGER(R_w),win.m_low,nPix,delta);
    /** normal vector defining intersection plane */
	SEXP R_plane = R_NilValue;
	PROTECT(R_plane = getAttrib(R_S, install("plane"))); ++nprotect;
	if(isNull(R_plane))
	 error(_("`plane` numeric vector must be provided as an attribute."));
	STGM::CVector3d n(REAL(R_plane));

	/** iterate over all intersection profiles */
	int type = 0;
	SEXP R_obj = R_NilValue;
	size_t num = (size_t) LENGTH(R_S);
	for(size_t k=0; k < num ; k++)
	{
	   PROTECT(R_obj = VECTOR_ELT(R_S,k));
	   type = INTEGER(VECTOR_ELT(R_obj,1))[0];   						// intersection type

	   /** Convert R object to the corresponding intersection type:
	    *  see 'convert_R_CylinderIntersections' for
	    *  indices of intersection profiles 'R_obj' */
	   if(type == STGM::ELLIPSE_2D)
	   {
		   /** CEllipse2 */
		   STGM::CEllipse2 ellipse2 = convert_C_Ellipse2(R_obj);
		   digitizer(ellipse2);

	   } else if(type == STGM::DISC) {
		   /** CCircle3 (disc)*/
		   STGM::CVector3d ctr(REAL(VECTOR_ELT(R_obj,2)));	 			// center of disc
		   STGM::CCircle3 disc3(ctr,REAL(VECTOR_ELT(R_obj,3))[0],n,k);	// radius
		   digitizer(disc3);

	   } else if(type == STGM::CAP) {
		   STGM::CVector3d ctr(REAL(VECTOR_ELT(R_obj,7)));   			// center of cap 'mPoint0'
		   STGM::CCircle3 disc3(ctr,REAL(VECTOR_ELT(R_obj,13))[0],n,k);
		   digitizer(disc3);

	   } else {
		   /** CEllipse3 including caps, segments, arcs as
		     * objects as part of the ellipse object in 3D */
		   STGM::CEllipse3 ellipse3 = convert_C_Ellipse3(R_obj, n);
		   digitizer(ellipse3);
	   }
	   UNPROTECT(1);
	}

	UNPROTECT(nprotect);
    return R_w;
}


namespace STGM
{

template<class T>
void CPoissonSystem<T>::simSystem(SEXP R_args, SEXP R_cond) {
	SEXP R_fname = R_NilValue, R_label = R_NilValue;
	PROTECT(R_fname = getListElement( R_cond, "rdist"));
	PROTECT(R_label = getListElement( R_cond, "label"));

	const char *label = "N";
    if(!isNull(R_label))
	  label = translateChar(asChar(R_label));
	int isPerfect = INTEGER(AS_INTEGER(getListElement( R_cond, "perfect" )))[0];

	if(TYPEOF(R_fname) != VECSXP)							/* user-defined simulation function */
	{
		SEXP R_call, R_rho;
		PROTECT(R_rho  = getListElement( R_cond,"rho" ) );
		PROTECT(R_call = getCall( R_fname, R_args, R_rho)  );
		const char *ftype = CHAR(STRING_ELT(R_fname,0));

		// simulate
		GetRNGstate();
		simJoint(R_call, R_rho, ftype, label);				/* non perfect */
		PutRNGstate();
		UNPROTECT(2);

	} else {

	    // distribution types
	    const char *ftype_size  = GET_NAME(R_fname,0);
	    const char *ftype_shape = GET_NAME(R_fname,1);
	    const char *ftype_dir   = GET_NAME(R_fname,2);

	    if( !std::strcmp( ftype_size, "rbinorm"))
	    {

	      // size parameters
	      SEXP R_tmp = R_NilValue;
	      PROTECT(R_tmp = VECTOR_ELT( R_args, 0));
	      if(isNull(R_tmp) || LENGTH(R_tmp) != 5)
	    	error(_("Number of arguments required for bivariate normal is invalid."));

	      double mx  = REAL(getListElement( R_tmp, "mx"))[0];
	      double my  = REAL(getListElement( R_tmp, "my"))[0];
	      double sdx = REAL(getListElement( R_tmp, "sdx"))[0];
	      double sdy = REAL(getListElement( R_tmp, "sdy"))[0];
	      double rho = REAL(getListElement( R_tmp, "rho"))[0];
	      UNPROTECT(1);

	      if(isPerfect) {

	    	 /// size
	    	 rbinorm_exact_t rdist(mx,my,sdx,sdy,rho,m_box,ftype_size);

	    	 if(PL>0) {
	    		 Rprintf("\n");
	    		 Rprintf("Size/Shape:  mx=%f, sdx=%f, my=%f, sdy=%f, rho=%f \n",rdist.mx,rdist.sdx,rdist.my,rdist.sdy,rdist.rho);
	    		 Rprintf("Cumulative sum of probabilities: %f, %f, %f, %f \n",rdist.p[0],rdist.p[1],rdist.p[2],rdist.p[3]);
	    	 }

	    	 GetRNGstate();

			 if(!std::strcmp(ftype_dir, "runifdir" ))
			 {
				 runidir_t rdirect;
				 simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

			 } else {

				// direction parameters
				PROTECT(R_tmp = getListElement( R_args, "orientation" ));
				if(isNull(R_tmp) || LENGTH(R_tmp) != 1)									/* only one parameter required */
				  error(_("Argument `orientation` must have length equal to one."));
				double kappa = REAL(VECTOR_ELT(R_tmp,0))[0];
				UNPROTECT(1);

				if(!std::strcmp( ftype_dir, "rbetaiso" )) {

					rbetaiso_t rdirect(m_mu,kappa);
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {

					rVonMisesFisher_t rdirect(m_mu,kappa);
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else if(!std::strcmp( ftype_dir, "const")) {

					constdir_t rdirect(kappa);
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else {
				   error(_("Direction distribution type is not supported."));
				}

			}

		    PutRNGstate();

	      }  else { /* end isPerfect */

	    	 /// size
	    	 isPerfect = 0;
			 rbinorm_t rdist(mx,my,sdx,sdy,rho,m_box,ftype_size);

			 if(PL>0) {
			   Rprintf("\n");
			   Rprintf("Size/Shape:  mx=%f, sdx=%f, my=%f, sdy=%f, rho=%f \n",rdist.mx,rdist.sdx,rdist.my,rdist.sdy,rdist.rho);
			 }

	    	 GetRNGstate();
			 if(!std::strcmp(ftype_dir, "runifdir" ))
			 {
				 runidir_t rdirect;
				 simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

			 } else {

				// direction parameters
				PROTECT(R_tmp = getListElement( R_args, "orientation" ));
				if(isNull(R_tmp) || LENGTH(R_tmp) != 1)
				  error(_("Argument `orientation` must have length equal to one."));
				double kappa = REAL(VECTOR_ELT(R_tmp ,0))[0];
				UNPROTECT(1);

				if(!std::strcmp( ftype_dir, "rbetaiso" )) {

					rbetaiso_t rdirect(m_mu,kappa);
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {

					rVonMisesFisher_t rdirect(m_mu,kappa);
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else {
					error(_("Direction distribution type is not supported."));
			 	}

			}

			PutRNGstate();
	      }

	    } else {

	    	/*  -Univariate distributions for size and shape,
	    	 *  -orientation distribution independent of size/shape
	    	 *  -always non perfect simulation, even for `rlnorm`
	    	 */

	    	/// size
	    	SEXP R_tmp = R_NilValue;
	    	PROTECT(R_tmp = getListElement( R_args, "size" ));
	       	if(isNull(R_tmp) || LENGTH(R_tmp) == 0)
	    	 error(_("Parameters for `size` cannot have length zero."));
	    	double p1 = REAL(VECTOR_ELT( R_tmp ,0))[0];
	    	double p2 = (LENGTH(R_tmp) > 1 ? REAL(VECTOR_ELT( R_tmp,1))[0] : 0.0);	  /* could be constant then only one parameter */
	    	UNPROTECT(1);

	    	// shape
	    	PROTECT(R_tmp = getListElement( R_args, "shape" ));
			if(isNull(R_tmp) || LENGTH(R_tmp) == 0)
			 error(_("Parameters for `shape` cannot have length zero."));
			double s1 = REAL(VECTOR_ELT( R_tmp ,0))[0];
			double s2 = (LENGTH(R_tmp) > 1 ? REAL(VECTOR_ELT( R_tmp,1))[0] : 0.0);
			UNPROTECT(1);

			if( isPerfect && !std::strcmp( ftype_size, "rlnorm")) {

				rlnorm_exact_t rsize(p1,p2,m_box,ftype_size);
				rndGen_t rshape(s1,s2, ftype_shape);
				rndUnivar_t<rlnorm_exact_t> rdist(rsize,rshape);

				if(PL>0) {
				  Rprintf("\n");
				  Rprintf("Size:  %f, %f \n",rsize.mx,rsize.sdx);
				  Rprintf("Shape: %f, %f \n",rshape.p1,rshape.p2);
				}

				GetRNGstate();
				if(!std::strcmp(ftype_dir, "runifdir" ))
				{
					runidir_t rdirect;
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else {

					// direction parameters
					PROTECT(R_tmp = getListElement( R_args, "orientation" ));
					if(isNull(R_tmp) || LENGTH(R_tmp) != 1)
					  error(_("Argument `orientation` must have length equal to one."));
					double kappa = REAL(VECTOR_ELT(R_tmp,0))[0];
					UNPROTECT(1);

					if(!std::strcmp( ftype_dir, "rbetaiso" )) {

						rbetaiso_t rdirect(m_mu,kappa);
						simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

					} else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {

						rVonMisesFisher_t rdirect(m_mu,kappa);
						simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

					} else {
					   error(_("Direction distribution type is not supported."));
					}
				}

				PutRNGstate();

			} else  {

				rndGen_t rsize(p1,p2,ftype_size,m_box.volume());
				rndGen_t rshape(s1,s2, ftype_shape);
				rndUnivar_t<rndGen_t> rdist(rsize,rshape);

				if(PL>0) {
					Rprintf("\n");
					Rprintf("Size:  %f, %f \n",rdist.rsize.p1,rdist.rsize.p2);
					Rprintf("Shape: %f, %f \n",rdist.rshape.p1,rdist.rshape.p2);
			    }

				GetRNGstate();
				if(!std::strcmp(ftype_dir, "runifdir" ))
				{
					runidir_t rdirect;
					simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

				} else {

					// direction parameters
					PROTECT(R_tmp = getListElement( R_args, "orientation" ));
					if(isNull(R_tmp) || LENGTH(R_tmp) != 1)
					  error(_("Argument `orientation` must have length equal to one."));
					double kappa = REAL(VECTOR_ELT(R_tmp,0))[0];
					UNPROTECT(1);

					if(!std::strcmp( ftype_dir, "rbetaiso" )) {

						rbetaiso_t rdirect(m_mu,kappa);
						simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

					} else if(!std::strcmp( ftype_dir, "rvMisesFisher")) {

						rVonMisesFisher_t rdirect(m_mu,kappa);
						simBivariate(rdist,rdirect,label,ftype_size,isPerfect);

					} else {
					   error(_("Direction distribution type is not supported."));
					}
				}

				PutRNGstate();
			}
	    }
	}

	if(PL>0){
	 Rprintf("Done. Simulated %d objects. \n", m_objects.size());
	}

	UNPROTECT(2);
}



/**
 * @brief simulation spheroid system but not perfect
 *        using user defined random joint distribution
 *
 */
template<>
void CPoissonSystem<CSpheroid>::simJoint(SEXP R_call, SEXP R_rho, const char* type, const char *label) {
     int nTry=0;
     double mu = m_box.volume();

     while(m_num==0 && nTry<MAX_ITER) {
       m_num = rpois(mu*m_lam);
       ++nTry;
     }
     m_objects.reserve(m_num);
     double *v=0, a=0, b=0, c=0, theta=0, phi=0;

     CVector3d u;
     SEXP Reval = R_NilValue;

     if(PL>0) {
    	 Rprintf("\n");
    	 Rprintf("Spheroid (joint) simulation with `%s` \n", type);
     }

     int info = 0;
     for (size_t niter=0; niter<m_num; niter++)
     {
         PROTECT(Reval = R_tryEval(R_call,R_rho,&info));
         if(info != 0){
           error(_("simJoint(): R 'try-error' in evaluation of user-defined distribution function."));
         }
         a = REAL(getListElement(Reval,"a"))[0];			// 1nd. semi-minor
         b = REAL(getListElement(Reval,"b"))[0];			// 		semi-major
         c = REAL(getListElement(Reval,"c"))[0];			// 2st. semi-minor
         theta = REAL(getListElement(Reval,"theta"))[0];
         phi = REAL(getListElement(Reval,"phi"))[0];

         v = REAL(getListElement(Reval,"u"));
         u[0] = v[0]; u[1] = v[1]; u[2] = v[2];

         CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
         		 	   	  runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
						  runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

         m_objects.push_back( CSpheroid(center,a,c,b,u,theta,phi,m_objects.size()+1,label) );
         UNPROTECT(1);
     }
}

/**
 * @brief Bivariate ellipsoid-size-shape distribution,
 *        the major semiaxis is logN distributed,
 *        and possibly shorter semiaxes of unequal lengths
 */
template<>
template<class T1, class DIR>
void CPoissonSystem<CSpheroid>::simBivariate(T1 &rdist, DIR &rdir, const char *label, const char *type, int perfect)
{
	/* get Poisson parameter */
	  int nTry=0;
	  double mu = rdist.mu;
	  while(m_num==0 && nTry<MAX_ITER) {
		m_num = rpois(mu*m_lam);
		++nTry;
	  }
	  m_objects.reserve(m_num);

	  if(PL>0)
	  {
		   Rprintf("\n");
		   Rprintf("Spheroid simulation with `%s` (perfect=%d): \n", type, perfect);
		   if(perfect)
			   Rprintf("Mean number: %f (exact simulation: %f) \n", m_lam, mu);
		   else
			   Rprintf("Mean number: %f (Box volume %f) \n", m_lam, mu);
		   Rprintf("Number of spheroids: %d \n", m_num);
		   Rprintf("Set label '%s'. \n",label);
		   Rprintf("\n\n");
	  }

      CVector3d u;
      double b=0,c=0,    							/*  shorter semi-minor lengths a,c and semi-major length b */
    		 s=1.0,phi=0,theta=0;

      if(perfect) {

    	  double r=0.0;
    	  for (size_t niter=0; niter<m_num; niter++)
    	  {
    		  rdist(s,b,c);							/* shape factor, semi-major, 1st. semi-minor */
    		  r=b;									/* copy because could be swapped with shorter axis */
    		  if(m_type==OBLATE)
    		    std::swap(b,c);

    		  rdir(u,theta,phi);					/* direction */

			  CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
			                   runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
			                   runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));

			  m_objects.push_back( CSpheroid(center,c,c,b,u,theta,phi,niter+1,label) );
    	  }

      } else {

    	  for (size_t niter=0; niter<m_num; niter++)
    	  {
    		  rdist(s,b,c);
			  if(m_type==OBLATE)
				std::swap(b,c);

			  rdir(u,theta,phi);					/* direction */
			  CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
			       		  	   runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
			       			   runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

			  m_objects.push_back( CSpheroid(center,c,c,b,u,theta,phi,niter+1,label) );

    	  }
      }

}



/**
 * @brief  User-defined joint simulation od cylinders
 */

template<>
void CPoissonSystem<CCylinder>::simJoint(SEXP R_call, SEXP R_rho, const char* type, const char *label) {
     int nTry=0;
     double mu = m_box.volume();
     while(m_num==0 && nTry<MAX_ITER) {
       m_num = rpois(mu*m_lam);
       ++nTry;
     }
     m_objects.reserve(m_num);
     double *v=0,h=0,theta=0,phi=0,r=0;

     CVector3d u;
     SEXP Reval = R_NilValue;

     if(PL>0) {
      	 Rprintf("\n");
      	 Rprintf("Cylinder (joint) simulation with `%s` \n", type);
     }

     int info = 0;
     for (size_t niter=0; niter<m_num; niter++)
     {
         Reval = R_tryEval(R_call,R_rho,&info);
         if(info != 0)
           error(_("simJoint(): R 'try-error' error in evaluation of user-defined distribution function."));

         h = REAL(getListElement(Reval,"h"))[0]; 				/* height of cylinder! */
         r = REAL(getListElement(Reval,"r"))[0];				/* radius */
         theta = REAL(getListElement(Reval,"theta"))[0];      /* polar angle */
         phi = REAL(getListElement(Reval,"phi"))[0];

         v = REAL(getListElement(Reval,"u"));
         u[0]=v[0]; u[1]=v[1];  u[2]=v[2];

         CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
              			  runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
              			  runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

         m_objects.push_back(CCylinder(center,u,h,r,theta,phi,niter+1,label) );
     }

}

template<>
template<class T1, class DIR>
void CPoissonSystem<CCylinder>::simBivariate(T1 &rdist, DIR &rdir, const char *label, const char *type, int perfect)
{
	/* get Poisson parameter */
	int nTry=0;
	double mu = rdist.mu;
	while(m_num==0 && nTry<MAX_ITER) {
		m_num = rpois(mu*m_lam);
		++nTry;
	}
	m_objects.reserve(m_num);

	if(PL>0)
	{
	   Rprintf("\n");
	   Rprintf("Cylinder simulation with `%s` (perfect=%d): \n", type, perfect);
	   if(perfect)
		   Rprintf("Mean number: %f (exact simulation: %f) \n", m_lam, mu);
	   else
		   Rprintf("Mean number: %f (Box volume %f) \n", m_lam, mu);
	   Rprintf("Number of cylinders: %d \n", m_num);
	   Rprintf("Set label '%s'. \n",label);
	   Rprintf("\n");
	}

    CVector3d u;
    double r=0, 					/* radius and half length */
    	   h=0,s=1,phi=0,theta=0;

    if(perfect)
    {
    	  double len2=0;			/* exact simulation radius */

    	  for (size_t niter=0; niter<m_num; niter++)
      	  {
      		  rdist(s,h,r);				/* h is full length of cylinder including caps */
			  len2 = 0.5*h;				/* half length is `radius` for exact simulation */
			  r /= 2.0;
			  h -= 2.0*r;  			    /* finally h is the height (excluding caps) */

			  // direction
			  rdir(u,theta,phi);

			  /* sample positions conditionally of radii distribution */
			  CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*len2)+(m_box.m_low[0]-len2),
							   runif(0.0,1.0)*(m_box.m_size[1]+2*len2)+(m_box.m_low[1]-len2),
							   runif(0.0,1.0)*(m_box.m_size[2]+2*len2)+(m_box.m_low[2]-len2));

			  m_objects.push_back(CCylinder(center,u,h,r,theta,phi,niter+1,label) );
      	  }

    } else {

     	  for (size_t niter=0; niter<m_num; niter++)
     	  {
     		 rdist(s,h,r);				/* h is full length of cylinder including caps */
     		 r /= 2.0;
     		 h -= 2.0*r;  			    /* finally h is the height (excluding caps) */

     		 // direction
     		 rdir(u,theta,phi);

			 /* sample positions conditionally of radii distribution */
     		 CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
     					   	  runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
     						  runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

			 m_objects.push_back(CCylinder(center,u,h,r,theta,phi,niter+1,label) );
     	  }
    }

}

void CPoissonSystem<CSphere>::simSystem(SEXP R_param, SEXP R_cond) {
	SEXP R_fname = R_NilValue, R_label = R_NilValue;
	PROTECT(R_fname = getListElement( R_cond, "rdist"));
	PROTECT(R_label = getListElement( R_cond, "label"));

	const char *label = "N";
	if(!isNull(R_label))
	  label = translateChar(asChar(R_label));
	int isPerfect = INTEGER(AS_INTEGER(getListElement( R_cond, "perfect" )))[0];

   if(TYPEOF(R_fname) != VECSXP)
   {
	   SEXP R_call, R_rho;
	   PROTECT(R_rho = getListElement( R_cond, "rho" ));
	   PROTECT(R_call = getCall(R_fname,R_param,R_rho));

	   //R_eval_t<double> reval(R_call,R_rho,m_box.volume());					/* specialization */
	   const char *ftype = CHAR(STRING_ELT(R_fname,0));

	   GetRNGstate();
	   simJoint(R_call, R_rho, ftype, label);

	   PutRNGstate();
	   UNPROTECT(2);

   } else {

	   SEXP R_tmp = R_NilValue;
	   PROTECT(R_tmp = getListElement( R_param, "size" ));
	   if(isNull(R_tmp) || LENGTH(R_tmp) == 0)
		   error(_("Parameters for `size` cannot have length zero."));

	   double p1 = REAL_ARG_LIST(R_tmp, 0);
	   double p2 = (LENGTH(R_tmp) > 1 ? REAL_ARG_LIST(R_tmp, 1) : 0.0);
	   UNPROTECT(1);

       /* radii distribution */
       const char *ftype = CHAR(STRING_ELT(VECTOR_ELT(R_fname,0), 0));

       GetRNGstate();
       if(isPerfect) {
    	   rlnorm_exact_t rdist(p1,p2,m_box,ftype);

    	   if(PL>0) {
    	    Rprintf("\t Size parameters: mx=%f, sdx=%f \n", rdist.mx, rdist.sdx);
    	 	Rprintf("\t Cumulative sum of probabilities: %f, %f, %f, %f ( mu=%f ) \n",rdist.p[0],rdist.p[1],rdist.p[2],rdist.p[3], rdist.mu);
    	   }

    	   simUnivar(rdist,label,ftype,isPerfect);

       } else {
    	   rndGen_t rdist(p1,p2,ftype,m_box.volume());
    	   if(PL>0) {
    		  Rprintf("\t Size parameters: mx=%f, sdx=%f \n", rdist.p1, rdist.p2);
    	   }

    	   simUnivar(rdist,label,ftype,isPerfect);
       }
       PutRNGstate();

   }

   if(PL>0){
	  Rprintf("Done. Simulated %d spheres.\n", m_objects.size());
   }

   UNPROTECT(2);
}


void CPoissonSystem<CSphere>::simJoint(SEXP R_call, SEXP R_rho, const char* type, const char *label)
{
	 int info=0, nTry=0;
	 double mu = m_box.volume();
	 while(m_num==0 && nTry<MAX_ITER) {
	   m_num = rpois(mu*m_lam);
	   ++nTry;
	 }
	 m_objects.reserve(m_num);

	 if(PL>0){
		 Rprintf("Spheres simulation by `%s`  \n", type);
		 Rprintf("Box volume: %f, lam: %f, number of spheres: %d \n", m_box.volume(), m_lam, m_num);
	 }

	 SEXP Reval = R_NilValue;
	 for (size_t niter=0; niter<m_num; niter++)
	 {
		 PROTECT(Reval = R_tryEval(R_call,R_rho,&info));
		 if(info != 0)
		  error(_("simJoint(): `try` error in user defined distribution function."));

		 CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
						  runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
						  runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

		 m_objects.push_back( CSphere(center,REAL(Reval)[0],m_objects.size()+1,label));
		 UNPROTECT(1);
	  }
}


template< typename F>
void CPoissonSystem<CSphere>::simUnivar(F &rsize, const char *label, const char *type, int perfect)
{
  /* get Poisson parameter */
  int nTry=0;
  double mu = rsize.mu;
  while(m_num==0 && nTry<MAX_ITER) {
	m_num = rpois(mu*m_lam);
	++nTry;
  }
  m_objects.reserve(m_num);

  if(PL>0)
  {
	 Rprintf("\n");
	 Rprintf("Spheres simulation with `%s` (perfect=%d):  \n", type, perfect);
	 if(perfect)
	    Rprintf("Mean number: %f (exact simulation: %f) \n", m_lam, mu);
	 else Rprintf("Mean number: %f (Box volume %f) \n", m_lam, mu);
	 Rprintf("Set label '%s'. \n",label);
	 Rprintf("\n");
  }

  double r = 0.0;

  if(perfect) {

	  Rprintf("perfect \n");
	  for (size_t niter=0; niter<m_num; niter++)
	  {
			r = rsize();

			CVector3d center(runif(0.0,1.0)*(m_box.m_size[0]+2*r)+(m_box.m_low[0]-r),
							 runif(0.0,1.0)*(m_box.m_size[1]+2*r)+(m_box.m_low[1]-r),
							 runif(0.0,1.0)*(m_box.m_size[2]+2*r)+(m_box.m_low[2]-r));


			m_objects.push_back( CSphere(center,r,m_objects.size()+1,label));
	  }


  } else {

	  for (size_t niter=0; niter<m_num; niter++)
	  {
			r = rsize();

			CVector3d center(runif(0.0,1.0)*(m_box.m_size[0])+(m_box.m_low[0]),
							 runif(0.0,1.0)*(m_box.m_size[1])+(m_box.m_low[1]),
							 runif(0.0,1.0)*(m_box.m_size[2])+(m_box.m_low[2]));

			m_objects.push_back( CSphere(center,r,m_objects.size()+1,label));
	  }

  }

}

/**
 * @brief Intersect only objects whose center is inside the window
 *
 */
template<class T>
void IntersectWithPlane(CPoissonSystem<T> &sp, typename Intersectors<T>::Type &intersected, SEXP R_cond)
{

  SEXP R_intern = R_NilValue;
  PROTECT(R_intern = getListElement(R_cond,"intern"));
  if(isNull(R_intern))
	error(_("`intern` attribute not found."));
  int intern = INTEGER(R_intern)[0];

  SEXP R_n = R_NilValue;
  PROTECT(R_n = getListElement(R_cond,"nsect"));
  if(isNull(R_n))
	error(_("`nsect` normal vector not found."));
  CVector3d n(REAL(R_n));
  UNPROTECT(2);

  /** set up plane */
  CPlane plane(n,asReal(getListElement(R_cond,"dz")));

  CBox3 &box = sp.box();
  typename CPoissonSystem<T>::objects_array_t &objects = sp.refObjects();

  /** assume left-down corner as the origin of the box */
  if(intern)
  {
	  int i=0, j=0;
	  plane.getPlaneIdx(i,j);

	  CWindow win(box.m_size[i],box.m_size[j]);
      for(size_t k=0; k < objects.size(); ++k) {
           Intersector<T> intersector( objects[k], plane, box.m_size);
           if(intersector.FindIntersection()) {
               if((intersector.getObject())->isInWindow(win))
                 intersected.push_back( intersector );
           }
       }

  } else {

	  for(size_t k=0; k < objects.size(); ++k) {
		 Intersector<T> intersector( objects[k], plane, box.m_size);
		 if(intersector.FindIntersection())
		   intersected.push_back( intersector );
	 }

  }

  if(PL>10){
    Rprintf("Plane normal to: [%f %f %f] \n", n[0],n[1],n[2]);
    Rprintf("Number of intersections: %d \n",intersected.size());
  }

}

} // end namespace


//--- converter functions ---

SEXP convert_R_Spheres(STGM::CPoissonSystem<STGM::CSphere> &sp)
{
  SEXP R_ret = R_NilValue;
  STGM::Spheres &spheres = sp.refObjects();
  size_t num = sp.size();

  if(PL>0){
  	Rprintf("Converting spheres ... \n");
  }

  if(PL == 10)
  {
	  /* return radii only */
	  PROTECT(R_ret = allocVector(REALSXP,num));
	  double *discs = REAL(R_ret);
	  for(size_t k=0; k<num; k++)
		discs[k] = spheres[k].r();

  } else {

	  PROTECT(R_ret = allocVector(VECSXP, num) );
	  SEXP R_tmp, R_center;
	  const char *nms[] = {"id", "center", "r", ""};

	  STGM::CBox3 &box = sp.box();
	  const STGM::LateralPlanes &planes = box.getLateralPlanes();

	  for(size_t k = 0; k < num; k++)
	  {
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
		STGM::CVector3d &center = sphere.center();
		SET_REAL_VECTOR(R_center,center);

		SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
		SET_VECTOR_ELT(R_tmp,1,R_center);
		SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

		setAttrib(R_tmp, install("label"), mkString(sphere.label()) );
		setAttrib(R_tmp, install("interior"), ScalarLogical(interior));
		setAttrib(R_tmp, install("area"), ScalarReal(sphere.projectionArea()));

		SET_VECTOR_ELT(R_ret,k,R_tmp);
		UNPROTECT(2);
	  }
  }

  SET_CLASS_NAME(R_ret,"spheres");
  UNPROTECT(1);
  return R_ret;
}


STGM::CEllipse2 convert_C_Ellipse2(SEXP R_E)
{
   STGM::CVector2d ctr(REAL(VECTOR_ELT(R_E,2)));
   STGM::CMatrix2d A(REAL(VECTOR_ELT(R_E,3)));

   return STGM::CEllipse2(A,ctr,INTEGER(VECTOR_ELT(R_E,0))[0]);

   /** Alternatively initialize 2D ellipse */

   /*
   double *ab = REAL(VECTOR_ELT(R_E,4));
   STGM::CVector2d minorA(REAL(VECTOR_ELT(R_E,5)));
   STGM::CVector2d majorA(REAL(VECTOR_ELT(R_E,6)));

   return STGM::CEllipse2(ctr,								// center
   		   	   	   	   	  A,								// matrix A
						  majorA,							// semi-major
						  minorA,							// semi-minor
						  ab[0],ab[1],						// a,b (semi-major, semi-minor length)
						  REAL(VECTOR_ELT(R_E,7))[0],		// phi
						  INTEGER(VECTOR_ELT(R_E,0))[0]);	// id


   */
}

STGM::Ellipses2 convert_C_Ellipses2(SEXP R_E) {
	STGM::Ellipses2 ellipses2;
	size_t num = (size_t) LENGTH(R_E);

	ellipses2.reserve(num);
	for(size_t i=0; i < num; i++)
	  ellipses2.push_back( convert_C_Ellipse2( VECTOR_ELT(R_E,i) ) );

	return ellipses2;
}


SEXP convert_R_Ellipse2(STGM::CEllipse2 &ellipse) {
  SEXP R_tmp = R_NilValue;
  SEXP R_center, R_minor, R_major, R_ab, R_A;
  const char *nms[] = {"id", "type", "center", "A", "ab", "minor", "major", "phi", "S", ""};

  PROTECT(R_tmp = mkNamed(VECSXP,nms));
  PROTECT(R_center = allocVector(REALSXP, 2));
  PROTECT(R_ab = allocVector(REALSXP, 2));
  PROTECT(R_A = allocMatrix(REALSXP,2,2));
  PROTECT(R_minor = allocVector(REALSXP, 2));
  PROTECT(R_major = allocVector(REALSXP, 2));

  STGM::CVector2d &center = ellipse.center();
  SET_REAL_VECTOR(R_center,center);

  STGM::CVector2d &minor = ellipse.minorAxis();
  SET_REAL_VECTOR(R_minor,minor);

  STGM::CVector2d &major = ellipse.majorAxis();
  SET_REAL_VECTOR(R_major,major);

  REAL(R_ab)[0] = ellipse.a();    // major semi-axis (for both prolate/oblate)
  REAL(R_ab)[1] = ellipse.b();	  // minor semi-axis (for both prolate/oblate)

  const STGM::CMatrix2d &A = ellipse.MatrixA();
  COPY_C2R_MATRIX(A,R_A,2);

  SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
  SET_VECTOR_ELT(R_tmp,1,ScalarInteger(STGM::ELLIPSE_2D));
  SET_VECTOR_ELT(R_tmp,2,R_center);
  SET_VECTOR_ELT(R_tmp,3,R_A);
  SET_VECTOR_ELT(R_tmp,4,R_ab);
  SET_VECTOR_ELT(R_tmp,5,R_minor);
  SET_VECTOR_ELT(R_tmp,6,R_major);

  /* return angle between [0,2pi] */
  SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.phi()));
  SET_VECTOR_ELT(R_tmp,8,ScalarReal(ellipse.b()/ellipse.a()));

  UNPROTECT(6);
  return(R_tmp);

}

/**
 * @brief Convert ellipses to R objects
 * 		  Need angle between [0,2pi] for
 * 		  plotting.
 */
SEXP convert_R_Ellipses(STGM::Intersectors<STGM::CSpheroid>::Type &objects, STGM::CBox3 &box) {
  SEXP R_ret = R_NilValue;
  PROTECT(R_ret = allocVector(VECSXP, objects.size()) );

  if( PL == 10)  									/* short version of ellipse list */
  {
	SEXP R_tmp = R_NilValue;
	const char *nms[] = {"A", "C", "S", "phi", ""};

	for(size_t k=0; k<objects.size(); ++k)
	{
		STGM::CEllipse2 &ellipse = objects[k].getEllipse();
		PROTECT(R_tmp = mkNamed(VECSXP, nms));
		SET_VECTOR_ELT(R_tmp,0,ScalarReal(ellipse.a()));                  // major semi-axis (for both prolate/oblate)
		SET_VECTOR_ELT(R_tmp,1,ScalarReal(ellipse.b()));                  // minor semi-axis (for both prolate/oblate)
		SET_VECTOR_ELT(R_tmp,2,ScalarReal(ellipse.b()/ellipse.a()));      // shape
		SET_VECTOR_ELT(R_tmp,3,ScalarReal(ellipse.phi()));		          // relative to ´x´ axis
		SET_VECTOR_ELT(R_ret,k,R_tmp);
		UNPROTECT(1);
	}

  } else {											/* full version of ellipse list */

	  SEXP R_tmp, R_center, R_minor, R_major, R_ab, R_A;
	  const char *nms[] = {"id", "type", "center", "A", "ab", "minor", "major", "phi", "S", ""};

	  for(size_t k=0; k<objects.size(); ++k)
	  {
		  STGM::CEllipse2 &ellipse = objects[k].getEllipse();

		  PROTECT(R_tmp = mkNamed(VECSXP,nms));
		  PROTECT(R_center = allocVector(REALSXP, 2));
		  PROTECT(R_ab = allocVector(REALSXP, 2));
		  PROTECT(R_A = allocMatrix(REALSXP,2,2));
		  PROTECT(R_minor = allocVector(REALSXP, 2));
		  PROTECT(R_major = allocVector(REALSXP, 2));

		  STGM::CVector2d &center = ellipse.center();
		  SET_REAL_VECTOR(R_center,center);

		  STGM::CVector2d &minor = ellipse.minorAxis();
		  SET_REAL_VECTOR(R_minor,minor);

		  STGM::CVector2d &major = ellipse.majorAxis();
		  SET_REAL_VECTOR(R_major,major);

		  REAL(R_ab)[0] = ellipse.a();    // major semi-axis (for both prolate/oblate)
		  REAL(R_ab)[1] = ellipse.b();	  // minor semi-axis (for both prolate/oblate)

		  const STGM::CMatrix2d &A = ellipse.MatrixA();
		  COPY_C2R_MATRIX(A,R_A,2);

		  SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
		  SET_VECTOR_ELT(R_tmp,1,ScalarInteger(STGM::ELLIPSE_2D));
		  SET_VECTOR_ELT(R_tmp,2,R_center);
		  SET_VECTOR_ELT(R_tmp,3,R_A);
		  SET_VECTOR_ELT(R_tmp,4,R_ab);
		  SET_VECTOR_ELT(R_tmp,5,R_minor);
		  SET_VECTOR_ELT(R_tmp,6,R_major);

		  /* return angle between [0,2pi] */
		  SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.phi()));
		  SET_VECTOR_ELT(R_tmp,8,ScalarReal(ellipse.b()/ellipse.a()));

		  SET_VECTOR_ELT(R_ret,k,R_tmp);
		  UNPROTECT(6);
	  }
  }

  // set intersecting plane normal vector
  SEXP R_plane;
  STGM::CPlane & plane = objects[0].getPlane();
  PROTECT(R_plane = allocVector(REALSXP,3));
  SET_REAL_VECTOR(R_plane, plane.n);
  setAttrib(R_ret,install("plane"),R_plane);

  SEXP R_win;
  PROTECT(R_win = allocVector(VECSXP,2));
  setWindow(R_win,box,plane);
  setAttrib(R_ret,install("win"),R_win);

  SET_CLASS_NAME(R_ret,"ellipses");
  UNPROTECT(3);
  return(R_ret);
}


STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  STGM::CVector3d ctr(REAL(VECTOR_ELT( R_spheroid, 1)));
  STGM::CVector3d   u(REAL(VECTOR_ELT( R_spheroid, 2)));

  double *acb = REAL(VECTOR_ELT(R_spheroid, 3));
  double *angles = REAL(VECTOR_ELT(R_spheroid, 4));

  return STGM::CSpheroid (ctr,acb[0],acb[1],acb[2],u,angles[0],angles[1],
			  INTEGER(VECTOR_ELT(R_spheroid, 0))[0],
			  translateChar(asChar(getAttrib(R_spheroid, install("label")))),
			  INTEGER(getAttrib(R_spheroid, install("interior")))[0]);
}


STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  STGM::Spheroids spheroids;
  size_t num = (size_t) LENGTH(R_spheroids);
  spheroids.reserve(num);

  for(size_t i=0; i < num; i++)  {
      spheroids.push_back( convert_C_Spheroid( VECTOR_ELT(R_spheroids,i) ) );
  }

  return spheroids;
}


STGM::CSphere convert_C_Sphere(SEXP R_sphere)
{
	STGM::CVector3d ctr(REAL(VECTOR_ELT( R_sphere, 1)));

	return STGM::CSphere(ctr,REAL(VECTOR_ELT(R_sphere, 2))[0],
			INTEGER(VECTOR_ELT(R_sphere, 0))[0],
			translateChar(asChar(getAttrib(R_sphere, install("label")))),
			INTEGER(getAttrib(R_sphere, install("interior")))[0]);
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres)
{
  STGM::Spheres spheres;
  size_t num = (size_t) LENGTH(R_spheres);
  spheres.reserve(num);

  for(size_t i=0; i < num; i++)
      spheres.push_back(convert_C_Sphere( VECTOR_ELT(R_spheres,i) ) );

  return spheres;
}


SEXP convert_R_Ellipsoids(STGM::CPoissonSystem<STGM::CSpheroid> &sp) {
  int nProtected=0, ncomps=6, dim=3;

  STGM::CBox3 &box = sp.box();
  STGM::Spheroids &spheroids = sp.refObjects();

  SEXP R_ret = R_NilValue;
  PROTECT(R_ret = allocVector(VECSXP, spheroids.size()) );
  ++nProtected;

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

    STGM::CVector3d &center = spheroid.center();
    SET_REAL_VECTOR(R_center,center);

    STGM::CVector3d &u = spheroid.u();
    SET_REAL_VECTOR(R_u,u);

    REAL(R_acb)[0]=spheroid.a();
    REAL(R_acb)[1]=spheroid.c();
    REAL(R_acb)[2]=spheroid.b();    // major semi-axis

    REAL(R_angles)[0]=spheroid.theta();
    REAL(R_angles)[1]=spheroid.phi();

    const STGM::CMatrix3d &M = spheroid.rotationMatrix();
    COPY_C2R_MATRIX(M,R_rotM,dim);

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

    SET_VECTOR_ELT(R_ret,k,R_tmp);
    UNPROTECT(7);
  }

  SET_CLASS_NAME(R_ret, sp.m_type_str);
  UNPROTECT(nProtected);
  return(R_ret);
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

  STGM::CVector3d &u = cyl.u();
  SET_REAL_VECTOR(R_u,u);
  STGM::CVector3d &center = cyl.center();
  SET_REAL_VECTOR(R_center,center);
  STGM::CVector3d &origin0 = cyl.origin0();
  SET_REAL_VECTOR(R_origin0,origin0);
  STGM::CVector3d &origin1 = cyl.origin1();
  SET_REAL_VECTOR(R_origin1,origin1);

  REAL(R_angles)[0]=cyl.theta();
  REAL(R_angles)[1]=cyl.phi();

  const STGM::CMatrix3d &M = cyl.rotationMatrix();
  COPY_C2R_MATRIX(M,R_rotM,dim);

  /*  check intersection with bounding box*/
  STGM::Intersector<STGM::CCylinder> intersector(cyl, box.m_size );
  Rboolean interior = (Rboolean) TRUE;
  for(size_t j=0; j<planes.size(); ++j) {
       if( intersector(planes[j])) {
         interior = (Rboolean) FALSE;
         break;
       }
  }
  /* projected area (Delamination) */
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
  SET_STRING_ELT(R_names, 4, mkChar("h"));
  SET_STRING_ELT(R_names, 5, mkChar("u"));
  SET_STRING_ELT(R_names, 6, mkChar("r"));
  SET_STRING_ELT(R_names, 7, mkChar("angles"));
  SET_STRING_ELT(R_names, 8, mkChar("rotM"));

  setAttrib(R_Cyl, R_NamesSymbol, R_names);
  setAttrib(R_Cyl, install("label"), mkString(cyl.label()) );
  setAttrib(R_Cyl, install("interior"), ScalarLogical(interior));
  setAttrib(R_Cyl, install("area"), ScalarReal(area));

  UNPROTECT(8);
  return R_Cyl;
}

STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  const char *label = "N";

  /* Ferrit - actually a sphere as a cylinder because of application of FBA
   * otherwise 'P' for particle and 'N' default (no label)
   * */
  SEXP R_label = R_NilValue;
  PROTECT(R_label = getAttrib(R_cyl, install("label")));
  if(!isNull(R_label))
    label = translateChar(asChar(R_label));
  else{
	error(_("Undefined attribute `label`."));
  }

  SEXP R_int = R_NilValue;
  PROTECT(R_int = getAttrib(R_cyl, install("interior")));
  if(!isNull(R_int))
    interior = INTEGER(R_int)[0];
  else {
    error(_("Undefined attribute `interior`."));
  }

  /** just copy from R */
  STGM::CVector3d ctr(REAL(VECTOR_ELT(R_cyl, 1))),
	   		            u(REAL(VECTOR_ELT(R_cyl, 5)));

  UNPROTECT(2);
  return STGM::CCylinder(ctr,u,
    		  	REAL(VECTOR_ELT(R_cyl,4))[0],
                REAL(VECTOR_ELT(R_cyl,6))[0],
				REAL(VECTOR_ELT(R_cyl,7))[0],
				REAL(VECTOR_ELT(R_cyl,7))[1],
                INTEGER(VECTOR_ELT(R_cyl, 0))[0],
				label, interior);
}


STGM::Cylinders convert_C_Cylinders(SEXP R_cyls)
{
  STGM::Cylinders cylinders;
  size_t num = (size_t) LENGTH(R_cyls);
  cylinders.reserve(num);

  for(size_t i=0; i<num; i++)
   cylinders.push_back( convert_C_Cylinder(VECTOR_ELT(R_cyls,i)));

  return cylinders;
}


void setWindow(SEXP R_win, STGM::CBox3 &box, STGM::CPlane &plane)
{
	// set intersection window
	int i=0, j=0;
	plane.getPlaneIdx(i,j);
	if(PL>0){
	 Rprintf("Getting plane indices: [%d %d ] \n", i, j);
	}
	SEXP R_winX, R_winY;
	PROTECT(R_winX = allocVector(REALSXP,2));
	PROTECT(R_winY = allocVector(REALSXP,2));

	REAL(R_winX)[0] = box.m_low[i];
	REAL(R_winX)[1] = box.m_up[i];
	REAL(R_winY)[0] = box.m_low[j];
	REAL(R_winY)[1] = box.m_up[j];

	SET_VECTOR_ELT(R_win,0,R_winX);
	SET_VECTOR_ELT(R_win,1,R_winY);
	UNPROTECT(2);
}



SEXP convert_R_Circles(STGM::Intersectors<STGM::CSphere>::Type & objects, STGM::CBox3 &box)
{
  SEXP R_ret = R_NilValue;
  size_t num = objects.size();
  if(PL>0){
	Rprintf("Converting %d discs. \n",num);
  }

  if(PL==10)
  {
	 /* return radii only */
	 PROTECT(R_ret = allocVector(REALSXP,num));
	 double *res = REAL(R_ret);
	 for(size_t k=0;k<num;k++)
	    res[k] = objects[k].getCircle().r();

  } else {

	  SEXP R_tmp, R_center;
	  const char *nms[] = {"id", "type", "center", "r", ""};
	  PROTECT(R_ret = allocVector(VECSXP, num) );

	  for(size_t k=0;k<num;k++)
	  {
		 STGM::CCircle3 &circle = objects[k].getCircle();
		 PROTECT(R_tmp = mkNamed(VECSXP, nms));

		 PROTECT(R_center = allocVector(REALSXP, 3));
		 STGM::CVector3d &center = circle.center();
		 SET_REAL_VECTOR(R_center, center);

		 /* type is needed here, though quite redundant, because of digitization */
		 SET_VECTOR_ELT(R_tmp,0,ScalarInteger(circle.Id()));
		 SET_VECTOR_ELT(R_tmp,1,ScalarInteger(STGM::DISC));
		 SET_VECTOR_ELT(R_tmp,2,R_center);
		 SET_VECTOR_ELT(R_tmp,3,ScalarReal(circle.r()));

		 SET_VECTOR_ELT(R_ret,k,R_tmp);
		 UNPROTECT(2);
	   }
  }

  SEXP R_win, R_plane;
  STGM::CPlane &plane = objects[0].getPlane();
  PROTECT(R_plane = allocVector(REALSXP,3));
  SET_REAL_VECTOR(R_plane, plane.n);
  setAttrib(R_ret,install("plane"),R_plane);

  PROTECT(R_win = allocVector(VECSXP,2));
  setWindow(R_win,box,plane);
  setAttrib(R_ret,install("win"),R_win);

  SET_CLASS_NAME(R_ret,"discs");
  UNPROTECT(3);
  return R_ret;
}


STGM::CEllipse3 convert_C_Ellipse3(SEXP R_E, STGM::CVector3d &n)
{
	STGM::CVector3d ctr(REAL(VECTOR_ELT(R_E,2)));
	STGM::CVector3d majorA(REAL(VECTOR_ELT(R_E,3)));
	STGM::CVector3d minorA(REAL(VECTOR_ELT(R_E,4)));

	/* midpoints of caps */
	STGM::CVector3d mPoint0(REAL(VECTOR_ELT(R_E,7)));
	STGM::CVector3d mPoint1(REAL(VECTOR_ELT(R_E,8)));

	double *ab    = REAL(VECTOR_ELT(R_E,9));
	double *psi   = REAL(VECTOR_ELT(R_E,12));
	double *rcaps = REAL(VECTOR_ELT(R_E,13));

	return STGM::CEllipse3(ctr,
					n,											// plane normal
					majorA, minorA,								// axes vectors
					mPoint0, mPoint1,							// midpoints of caps
					ab[0],ab[1],								// lengths
					REAL(VECTOR_ELT(R_E,10))[0],				// angle phi
					psi[0],psi[1],								// segment angles
					rcaps[0],rcaps[1],
					INTEGER(VECTOR_ELT(R_E,1))[0],				// intersection type
					INTEGER(VECTOR_ELT(R_E,14))[0]);			// `m_side` to set reference side `m_side0`.
}


SEXP convert_R_Cylinders( STGM::CPoissonSystem<STGM::CCylinder> &sp )
{
  STGM::Cylinders &cyls = sp.refObjects();
  size_t num = cyls.size();
  if(PL>0){
   Rprintf("Convert %d cylinders.\n", num);
  }

  SEXP R_ret = R_NilValue;
  PROTECT(R_ret = allocVector(VECSXP, num) );

  if(PL == 10)
  {
	    const char *nms[] = {"id", "h", "r", "angles", ""};
	    SEXP R_tmp, R_angles;

	    for(size_t k=0;k<num;++k)
	    {
	  	    PROTECT(R_tmp = mkNamed(VECSXP, nms));
	  	    PROTECT(R_angles = allocVector(REALSXP,2));

	  	    STGM::CCylinder &cyl = cyls[k];
	  	    REAL(R_angles)[0]=cyl.theta();
	  	    REAL(R_angles)[1]=cyl.phi();

	  	    SET_VECTOR_ELT(R_tmp,0,ScalarInteger( (int) cyl.Id() ));
	  	    SET_VECTOR_ELT(R_tmp,1,ScalarReal(cyl.h())) ;
	  	    SET_VECTOR_ELT(R_tmp,2,ScalarReal(cyl.r()));
	  	    SET_VECTOR_ELT(R_tmp,3,R_angles);

	  	    SET_VECTOR_ELT(R_ret,k,R_tmp);
	  	    UNPROTECT(2);
	    }

  } else {
	  // get lateral bounding planes
	  STGM::CBox3 &box = sp.box();
	  STGM::LateralPlanes &planes = box.getLateralPlanes();

	  for(size_t k=0;k<num;++k)
	 	 SET_VECTOR_ELT(R_ret,k,convert_R_Cylinder(cyls[k],planes,box));
  }

  SET_CLASS_NAME(R_ret,"cylinders");
  UNPROTECT(1);
  return(R_ret);
}

SEXP convert_R_CylinderIntersections(STGM::Intersectors<STGM::CCylinder>::Type &objects, STGM::CBox3 &box)
{
  int  nLoopProtected=0, ncomps=15, ncompsCircle=4;
  SEXP R_result, R_center, R_minor, R_major, R_ipt0, R_ipt1,
       R_mPoint0, R_mPoint1, R_height, R_ab, R_rcaps, R_psi,
	   R_obj = R_NilValue, R_names = R_NilValue;

  int type  = 0;
  size_t num = objects.size();
  PROTECT(R_result = allocVector(VECSXP, num) );

  for(size_t i=0; i<num; ++i)
  {
	  type = objects[i].getType();

      if(type == STGM::ELLIPSE ||
         type == STGM::ELLIPSE_ARC ||
         type == STGM::ELLIPSE_SEGMENT ||
         type == STGM::CAP) {

    	  nLoopProtected=0;
          PROTECT(R_obj = allocVector(VECSXP, ncomps) ); ++nLoopProtected;
          PROTECT(R_names = allocVector(STRSXP, ncomps));++nLoopProtected;

          PROTECT(R_center  = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_minor   = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_major   = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_ipt0    = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_ipt1    = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_mPoint0 = allocVector(REALSXP, 3) ); ++nLoopProtected;
          PROTECT(R_mPoint1 = allocVector(REALSXP, 3) ); ++nLoopProtected;

          PROTECT(R_ab = allocVector(REALSXP, 2) );      ++nLoopProtected;
          PROTECT(R_psi = allocVector(REALSXP, 2) );     ++nLoopProtected;
          PROTECT(R_rcaps = allocVector(REALSXP, 2) );   ++nLoopProtected;
          PROTECT(R_height = allocVector(REALSXP, 2) );  ++nLoopProtected;

          STGM::CEllipse3 &ellipse = objects[i].getEllipse();
          STGM::CVector3d &center = ellipse.center();
          SET_REAL_VECTOR(R_center,center);

          STGM::CVector3d &majorAxis = ellipse.majorAxis();
          SET_REAL_VECTOR(R_major, majorAxis);

          STGM::CVector3d &minorAxis = ellipse.minorAxis();
          SET_REAL_VECTOR(R_minor, minorAxis);

          STGM::CVector3d &ipt0 = objects[i].ipt0;
          SET_REAL_VECTOR(R_ipt0, ipt0);

          STGM::CVector3d &ipt1 = objects[i].ipt1;
          SET_REAL_VECTOR(R_ipt1, ipt1);

          REAL(R_ab)[0] = ellipse.a();
          REAL(R_ab)[1] = ellipse.b();
          REAL(R_psi)[0]= ellipse.psi()[0];
          REAL(R_psi)[1]= ellipse.psi()[1];

          /* no conversion to [0,pi/2] because of plotting */
          double phi = objects[i].getCylinder().phi();

          STGM::CVector3d &mPoint0 = objects[i].getCircle1().center();
          SET_REAL_VECTOR(R_mPoint0,mPoint0);

          STGM::CVector3d &mPoint1 = objects[i].getCircle2().center();
          SET_REAL_VECTOR(R_mPoint1,mPoint1);

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

      } else if(type == STGM::DISC ) {

          PROTECT(R_obj = allocVector(VECSXP, ncompsCircle) );   ++nLoopProtected;
          PROTECT(R_names = allocVector(STRSXP, ncompsCircle));  ++nLoopProtected;
          PROTECT(R_mPoint0 = allocVector(REALSXP, 3) );         ++nLoopProtected;

          /* circle1 stores the circle as intersection of cylinder */
          STGM::CVector3d &mPoint0 = objects[i].getCircle1().center();
          SET_REAL_VECTOR(R_mPoint0,mPoint0);

          SET_VECTOR_ELT(R_obj,2, R_mPoint0 );
          SET_VECTOR_ELT(R_obj,3,ScalarReal( objects[i].getCircle1().r()));
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
  }

  SEXP R_plane;
  STGM::CPlane & plane = objects[0].getPlane();
  PROTECT(R_plane = allocVector(REALSXP,3));
  SET_REAL_VECTOR(R_plane, plane.n);
  setAttrib(R_result,install("plane"),R_plane);

  SEXP R_win;
  PROTECT(R_win = allocVector(VECSXP,2));
  setWindow(R_win,box,plane);
  setAttrib(R_result,install("win"),R_win);

  SET_CLASS_NAME(R_result,"cylsects");
  UNPROTECT(3);
  return R_result;
}

