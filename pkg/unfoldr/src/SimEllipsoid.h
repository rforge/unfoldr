/**
 * SimEllipsoid.h
 *
 *  Created on: 10.01.2014
 *      Author: franke
 */

#ifndef SIMELLIPSOID_H_
#define SIMELLIPSOID_H_

#define MAX_ITER 100

#include "Utils.h"
#include "Intersector.h"

#ifdef __cplusplus
extern "C" {
#endif
 //
 SEXP EllipsoidSystem(SEXP R_param, SEXP R_cond);

 SEXP SimulateSpheroidsAndIntersect(SEXP R_param, SEXP R_cond, SEXP R_n);

 SEXP IntersectSpheroidSystem(SEXP R_S, SEXP R_n, SEXP R_z, SEXP R_intern, SEXP R_env, SEXP R_pl);

 SEXP DigitizeEllipseIntersections(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_delta);

 SEXP UpdateIntersections(SEXP RS, SEXP R_box);

#ifdef __cplusplus
}
#endif

namespace STGM {


/**
 *  Ellipsoid or Spheroid system
 */
class CSpheroidSystem
{
 public:
  CSpheroid::spheroid_type m_stype;

  CSpheroidSystem(CBox3 &box, double lam, CVector3d &mu, CSpheroid::spheroid_type stype, int perfect = 1)
    :  m_stype(stype), m_box(box), m_lam(lam), m_maxR(0), m_mu(mu),  num(0), m_perfect(perfect)
  {
  };

  ~CSpheroidSystem() {};


  void simSpheroidSystem(SEXP R_param, SEXP R_cond);

  void simUnivar(SEXP R_args, rdist2_t rsize, rdist2_t rshape, STGM::CSpheroid::direction_type dtype, const char *label);

  void simJoint(SEXP R_call, SEXP R_rho, const char *label);

  void simBivariate(SEXP R_args, rdist2_t rshape,STGM::CSpheroid::direction_type dtype, const char *label, int perfect);

  STGM::Spheroids &refObjects()  { return m_spheroids; }
  const STGM::Spheroids &refObjects() const { return m_spheroids; }

  inline size_t size()  { return m_spheroids.size(); }
  inline double maxR()  { return m_maxR; }
  inline double lam()   { return m_lam; }

  const STGM::CBox3 &box() const { return m_box; }
  STGM::CBox3 &box()  { return m_box; }

  CVector3d & u() { return m_mu; }
  const CVector3d u() const { return m_mu; }

  void IntersectWithPlane(STGM::Intersectors<STGM::CSpheroid>::Type &objects, STGM::CPlane &plane, int intern);

  Rboolean isPerfect() { return (m_perfect > 0 ? TRUE : FALSE); }

 private:
  CBox3 m_box;

  double m_lam, m_maxR;
  CVector3d m_mu;

  Spheroids m_spheroids;
  size_t num;
  int m_perfect;
};



} /* namespace STGM */

#endif /* SIMELLIPSOID_H_ */
