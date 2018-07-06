/*
 * simSphere.h
 *
 *  Created on: 07.05.2015
 *      Author: franke
 */

#ifndef SRC_SIMSPHERE_H_
#define SRC_SIMSPHERE_H_

#include "Utils.h"
#include "Rheaders.h"
#include "Intersector.h"

#ifdef __cplusplus
extern "C" {
#endif

  SEXP SphereSystem(SEXP R_param, SEXP R_cond);

  SEXP SimulateSpheresAndIntersect(SEXP R_param, SEXP R_cond);

  SEXP IntersectSphereSystem(SEXP R_var, SEXP R_n, SEXP R_z, SEXP R_intern, SEXP R_env, SEXP R_pl);

  SEXP DigitizeDiscs(SEXP R_S, SEXP R_cond, SEXP R_delta);

#ifdef __cplusplus
}
#endif


namespace STGM {

class CBoolSphereSystem {
public:

  CBoolSphereSystem(CBox3 &box, double lam, int perfect) :
    m_box(box), m_lam(lam), num(0), m_perfect(perfect)
  {}

  ~CBoolSphereSystem() {};

  size_t getNumSpheres() const { return num; }
  STGM::Spheres &refObjects()  { return m_spheres; }
  const STGM::Spheres &refObjects() const { return m_spheres; }

  void simSphereSys(SEXP R_param, SEXP R_cond);

  void IntersectWithPlane(STGM::Intersectors<STGM::CSphere>::Type &objects, STGM::CPlane &plane, int intern);

  template<typename F> void simSpheres(F f, const char *label);
  void simSpheresPerfect(double mx, double sdx, const char *label, int perfect);

  STGM::CBox3 & box() { return m_box; }

private:
  CBox3 m_box;
  double m_lam;
  size_t num;
  int m_perfect;
  Spheres m_spheres;

};

}


#endif /* SRC_SIMSPHERE_H_ */
