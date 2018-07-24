/**
 * SimEllipsoid.h
 *
 *  Created on: 10.01.2014
 *      Author: franke
 */

#ifndef SIMPOISSON_H_
#define SIMPOISSON_H_

#define MAX_ITER 100

#include "Utils.h"
#include "directions.h"
#include "Intersector.h"

#ifdef __cplusplus
extern "C" {
#endif

 SEXP PoissonSystem(SEXP R_param, SEXP R_cond);

 SEXP DigitizeProfiles(SEXP R_var, SEXP R_delta, SEXP R_win, SEXP R_env);

 SEXP IntersectPoissonSystem(SEXP R_var, SEXP R_cond, SEXP R_env);

 SEXP UpdateIntersections(SEXP R_S, SEXP R_env);


#ifdef __cplusplus
}
#endif



namespace STGM {


/**
 *  A class for general Poisson system
 */

template<class T>
class CPoissonSystem
{

 typedef T poisson_t;

 //typedef typename ConverterFunction<T>::Type converter_t;

 public:

  typedef std::vector<poisson_t> objects_array_t;

  const char * m_type_str;
  typedef enum { PROLATE = 0, OBLATE = 1, SPHERE = 2, CYLINDER = 3 } grain_type;
  typedef enum { UNIFORM_D = 0, BETAISOTROP_D = 1, MISES_D = 2} direction_type;

  CPoissonSystem(CBox3 &box, double lam, CVector3d &mu, const char * type, int perfect = 1)
    :  m_type_str(type), m_box(box), m_lam(lam), m_mu(mu), m_num(0), m_perfect(perfect)
  {
	  if( !std::strcmp(type, "prolate" ))
	   { m_type = PROLATE;}
	  else if(!std::strcmp(type, "oblate" ))
	   { m_type = OBLATE; }
	  else if(!std::strcmp(type, "cylinder" ))
	   { m_type = CYLINDER; }
	  // defaults
	  m_dirtype = UNIFORM_D;

	  // init converter
	  //m_converter = ConverterFunction<T>();

  };

  ~CPoissonSystem() {};

  void simSystem(SEXP R_param, SEXP R_cond);

  void simJoint(SEXP R_call, SEXP R_rho, const char *label);

  template<class T1, class DIR>
  void simBivariate(T1 &rdist, DIR &rdir, const char *label, const char *type, int perfect);

  objects_array_t &refObjects()  { return m_objects; }
  const objects_array_t &refObjects() const { return m_objects; }

  inline size_t size() { return m_objects.size(); }
  inline double lam() { return m_lam; }

  const CBox3 &box() const { return m_box; }
  CBox3 &box()  { return m_box; }

  CVector3d & u() { return m_mu; }
  const CVector3d u() const { return m_mu; }

  Rboolean isPerfect() { return (m_perfect > 0 ? TRUE : FALSE); }

  //ConverterFunction<T> &getConverter() { return m_converter; }

 private:

  CBox3 m_box;

  double m_lam;
  CVector3d m_mu;

  objects_array_t m_objects;
  size_t m_num;
  int m_perfect;
  grain_type m_type;
  direction_type m_dirtype;

  //ConverterFunction<T> m_converter;

};


template<>
class CPoissonSystem<CSphere>
{
public:

  typedef std::vector<CSphere> objects_array_t;

  const char * m_type_str;
  typedef enum { SPHERE = 0 } grain_type;

  CPoissonSystem(CBox3 &box, double lam, CVector3d &mu, const char * type,  int perfect) :
	  m_type_str(type), m_box(box), m_lam(lam), m_mu(mu), m_num(0), m_perfect(perfect), m_type(SPHERE)
  {
	  // init converter
	  //m_converter = ConverterFunction<CSphere>();
  }

  ~CPoissonSystem() {};

  size_t size() const { return m_objects.size(); }

  Spheres &refObjects()  { return m_objects; }
  const Spheres &refObjects() const { return m_objects; }

  void simSystem(SEXP R_param, SEXP R_cond);

  template<typename F>
  void simUnivar(F &size, const char *label, const char *type, int perfect);

  void IntersectWithPlane(Intersectors<CSphere>::Type &objects, CPlane &plane, int intern);

  CBox3 &box()  { return m_box; }
  const CBox3 &box() const { return m_box; }

  Rboolean isPerfect() { return (m_perfect > 0 ? TRUE : FALSE); }

  //ConverterFunction<CSphere> &getConverter() { return m_converter; }

private:
  CBox3 m_box;
  double m_lam;
  CVector3d m_mu;
  size_t m_num;

  int m_perfect;
  Spheres m_objects;
  grain_type m_type;

  //ConverterFunction<CSphere> m_converter;
};


struct runidir_t
{
   void operator()(CVector3d &u, double &theta, double &phi) { runidir(u.ptr(),theta,phi); }
};

struct rbetaiso_t
{
 CVector3d mu;
 double kappa;

 rbetaiso_t(CVector3d _mu, double _kappa) : mu(_mu), kappa(_kappa)
 {};

 void operator()(CVector3d &u, double &theta, double &phi)
 {
	 if(kappa < 1e-8) {
	    u = (runif(0.0,1.0)<0.5) ? mu : -mu;
     } else {
        rOhserSchladitz(u.ptr(),mu.ptr(),kappa,theta,phi);
	 }
 }

};

struct rVonMisesFisher_t
{
 CVector3d mu;
 double kappa;

 rVonMisesFisher_t(CVector3d _mu, double _kappa) : mu(_mu), kappa(_kappa)
 {};

 void operator()(CVector3d &u, double &theta, double &phi)
 {
	 if(kappa<1e-8)
	   runidir(u.ptr(),theta,phi);
     else
	   rVonMisesFisher(u.ptr(),mu.ptr(),kappa,theta,phi);
 }

};



struct rndGen_t {
  double p1, p2;
  rdist2_t rdist;
  double mu;

  rndGen_t(double _p1, double _p2, double _mu, const char* ftype)
    : p1(_p1), p2(_p2), mu(_mu)
  {
    if (!std::strcmp(ftype, "rbeta" )) {
        rdist = &rbeta;
    } else if(!std::strcmp(ftype, "rlnorm")) {
    	rdist= &rlnorm;
    } else if(!std::strcmp(ftype, "rgamma")) {
    	rdist = &rgamma;
    } else if(!std::strcmp(ftype, "runif" )) {
    	rdist = &runif;
    } else if(!std::strcmp(ftype, "const" )) {
    	rdist = &rconst;
    } else {
        error("Undefined `size` distribution function.");
    }

  };

  double operator()() { return rdist(p1,p2); }

};

struct rndSizeShape_t {
  rndGen_t rsize, rshape;
  double mu;

  rndSizeShape_t(double p1, double p2, double s1, double s2,
		  	  	   double _mu, const char *size, const char *shape)

    :  rsize(rndGen_t(p1,p2,_mu,size)), rshape(rndGen_t(s1,s2,_mu,shape)), mu(_mu)
  {
  };

  void operator()(double &s, double &b, double &c) {
	 b = rsize();
	 c = rshape();
	 s = c/b;
  }

};


/* only sphere radius  */
struct rlnorm_exact_t {
	double mx, sdx;
	const char *size;
	double mu, sdx2, p[4];


	rlnorm_exact_t(double _mx, double _sdx, CBox3 &box, const char *_size) :
		mx(_mx), sdx(_sdx), size(_size), mu(0)
	{
		sdx2 = SQR(sdx);
		// calculates `mu`
		cum_prob_k(mx,sdx2,box.m_up[0],box.m_up[1],box.m_up[2],p,&mu);
	}

	double operator()() {
		int k = sample_k(p);
		return rlnorm(mx+k*sdx2,sdx);
	}

};

/* bivariate normal using exact simulation */
struct rbinorm_exact_t {
	double mx,my,sdx,sdy,rho;
	double mu, sdx2, p[4], x, y;
	const char *size;

	rbinorm_exact_t(double _mx, double _my, double _sdx, double _sdy, double _rho,
						CBox3 &box, const char *_size)

	  :	mx(_mx), my(_my), sdx(_sdx), sdy(_sdy), rho(_rho), mu(0), x(0), y(0), size(_size)
	{
		sdx2 = SQR(sdx);
		cum_prob_k(mx,sdx2,box.m_up[0],box.m_up[1],box.m_up[2],p,&mu);
	}

	void operator()(double &s, double &b, double &c) {
		rbinorm_exact(p,mx,sdx,my,sdy,rho,x,y);
		s=1.0/(1.0+std::exp(-y));
		b=std::exp(x); 						/* a = r for exact simulation*/
		c=b*s;
	}

};

/* bivariate normal non exact simulation */
struct rbinorm_t {
	double mx,my,sdx,sdy,rho;
	double x, y;
	const char *size;
	double mu;

	rbinorm_t(double _mx, double _my, double _sdx, double _sdy, double _rho, double _mu, const char *_size)
		: mx(_mx), my(_my), sdx(_sdx), sdy(_sdy), rho(_rho),
		  x(0), y(0), size(_size), mu(_mu)
	{
	};

	void operator()(double &s, double &b, double &c) {
		rbinorm(mx,sdx,my,sdy,rho,x,y);
		s=1.0/(1.0+std::exp(-y));
		b=std::exp(x); 						/* a = r for exact simulation*/
		c=b*s;
	}

};




} /* namespace STGM */


/*
template<class T>
struct ConverterFunction { ConverterFunction(); };

template<>
struct ConverterFunction<CSpheroid> {
  typedef Spheroids Type;
  typedef Intersectors<CSpheroid>::Type intersector_t;

  //intersector_t get2DObjects(SEXP R_sp) { };

  SEXP set2DObjects(intersector_t &objects) { return R_NilValue; };

  Type get3DObjects(SEXP Rs) { return convert_C_Spheroids(Rs); };

  SEXP set3DObjects(CPoissonSystem<CSpheroid> &sp, int pl){ return convert_R_Ellipsoids(sp); };

};

template<>
struct ConverterFunction<CCylinder> {
  typedef CCylinder Type;
  typedef Intersectors<CCylinder>::Type intersector_t;


    //intersector_t get2DObjects(SEXP R_sp) { };

    SEXP set2DObjects(intersector_t &objects) { return R_NilValue; };

    Type get3DObjects(SEXP Rs) { return convert_C_Cylinders(Rs); };

    SEXP set3DObjects(CPoissonSystem<CCylinder> &sp) { return convert_R_Cylinders(sp); };
};


template<>
struct ConverterFunction<CSphere> {
  typedef CSphere Type;
  typedef Intersectors<CSphere>::Type intersector_t;


    //intersector_t get2DObjects(SEXP R_sp)  { };

    SEXP set2DObjects(intersector_t &objects) { return R_NilValue; };

    Type get3DObjects(SEXP Rs)   {	  return convert_C_Spheres(Rs);  };

    SEXP set3DObjects(CPoissonSystem<CSphere> &sp, int pl) { return convert_R_Spheres(sp); };
};

*/

#endif /* SIMPOISSON_H_ */
