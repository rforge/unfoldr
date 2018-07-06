/**
 * @file SimCylinder.h
 *
 *  @date: 09.05.2016
 *  @author: M. Baaske
 */

#ifndef SRC_SIMCYLINDER_H_
#define SRC_SIMCYLINDER_H_

#define MAX_ITER 100

#include "Utils.h"
#include "Intersector.h"

#ifdef __cplusplus
extern "C" {
#endif

  SEXP CylinderSystem(SEXP R_param, SEXP R_cond);

  SEXP IntersectCylinderSystem(SEXP R_var, SEXP R_n, SEXP R_dz, SEXP R_intern, SEXP R_env, SEXP R_pl);

  SEXP SimulateCylindersAndIntersect(SEXP R_param, SEXP R_cond);

#ifdef __cplusplus
}
#endif


namespace STGM {

class CCylinderSystem {
 public:

  CCylinderSystem(CBox3 &box, double lam, CVector3d &mu, int perfect = 1 )
      : m_box(box), m_lam(lam), m_maxR(0), m_mu(mu), num(0), m_perfect(perfect)
  {
    //box.ConstructBoundingPlanes();
  }

  void simCylinderSystem(SEXP R_param, SEXP R_cond);

  void simUnivar(SEXP R_args, rdist2_t rsize, rdist2_t rshape, STGM::CCylinder::direction_type dtype, const char *label);

  void simJoint(SEXP R_call, SEXP R_rho, const char *label);

  void simBivariate(SEXP R_args, STGM::CCylinder::direction_type dtype, const char *label, int perfect);

  inline size_t size()  { return m_cylinders.size(); }
  inline double maxR()  { return m_maxR; }

  Cylinders &refObjects()  { return m_cylinders; }
  const Cylinders &refObjects() const { return m_cylinders; }

  STGM::CBox3 &box()  { return m_box; }
  const STGM::CBox3 &box() const { return m_box; }

  void IntersectWithPlane(STGM::Intersectors<STGM::CCylinder>::Type &objects, STGM::CPlane &plane, int intern);
  Rboolean isPerfect() { return (m_perfect > 0 ? TRUE : FALSE); }

  private:
    CBox3 m_box;
    double m_lam, m_maxR;
    CVector3d m_mu;
    size_t num;

    Cylinders m_cylinders;
    int m_perfect;

};

} // STGM

#endif /* SRC_SIMCYLINDER_H_ */
