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

  void simUnivar(SEXP R_args, rdist2_t rsize, rdist2_t rshape,
		  	  	  direction_type dtype, const char *label);

  void simJoint(SEXP R_call, SEXP R_rho, const char *label);

  void simBivariate(SEXP R_args, rdist2_t rshape, direction_type dtype,
		             const char *label, int perfect);

  objects_array_t &refObjects()  { return m_objects; }
  const objects_array_t &refObjects() const { return m_objects; }

  inline size_t size()  { return m_objects.size(); }
  inline double lam()   { return m_lam; }

  const STGM::CBox3 &box() const { return m_box; }
  STGM::CBox3 &box()  { return m_box; }

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
class CPoissonSystem<STGM::CSphere>
{
public:

  typedef std::vector<STGM::CSphere> objects_array_t;

  const char * m_type_str;
  typedef enum { SPHERE = 0 } grain_type;

  CPoissonSystem(CBox3 &box, double lam,CVector3d &mu, const char * type,  int perfect) :
	  m_type_str(type), m_box(box), m_lam(lam), m_mu(mu), m_num(0), m_perfect(perfect), m_type(SPHERE)
  {
	  // init converter
	  //m_converter = ConverterFunction<CSphere>();
  }

  ~CPoissonSystem() {};

  size_t size() const { return m_objects.size(); }

  STGM::Spheres &refObjects()  { return m_objects; }
  const STGM::Spheres &refObjects() const { return m_objects; }

  void simSystem(SEXP R_param, SEXP R_cond);

  template<typename F> void simUnivar(F f, const char *label);
  void simSpheresPerfect(double mx, double sdx, const char *label, int perfect);

  void IntersectWithPlane(STGM::Intersectors<STGM::CSphere>::Type &objects, STGM::CPlane &plane, int intern);

  STGM::CBox3 &box()  { return m_box; }
  const STGM::CBox3 &box() const { return m_box; }

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
