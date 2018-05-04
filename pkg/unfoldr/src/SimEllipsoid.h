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
 SEXP IntersectSpheroidSystem(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_intern, SEXP R_pl);
 SEXP DigitizeEllipseIntersections(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_delta);
 SEXP GetEllipsoidSystem(SEXP ext);
 SEXP GetMaxRadius(SEXP ext);
 SEXP SetupSpheroidSystem(SEXP R_vname, SEXP R_env, SEXP R_param, SEXP R_cond);
 SEXP CopySpheroidSystem(SEXP R_spheroids, SEXP R_env, SEXP R_param, SEXP R_cond);
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

  void simSpheroidSys(R_Calldata d);
  void simConstSpheroidSys(R_Calldata d);
  void simSysJoint(R_Calldata d);
  void simBivariate(R_Calldata d);
  void simBivariate2(R_Calldata d);  									/* non equal shorter semiaxes */

  STGM::Spheroids &refObjects()  { return m_spheroids; }
  const STGM::Spheroids &refObjects() const { return m_spheroids; }

  inline size_t size()  { return m_spheroids.size(); }
  inline double maxR()  { return m_maxR; }

  const STGM::CBox3 &box() const { return m_box; }
  STGM::CBox3 &box()  { return m_box; }

  void IntersectWithPlane(STGM::Intersectors<STGM::CSpheroid>::Type &objects, STGM::CPlane &plane, int intern);

  Rboolean isPerfect() { return (m_perfect > 0 ? TRUE : FALSE); }

 private:
  CBox3 m_box;

  double m_lam,m_maxR;
  CVector3d m_mu;

  Spheroids m_spheroids;
  size_t num;
  int m_perfect;
};



} /* namespace STGM */

#endif /* SIMELLIPSOID_H_ */
