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

  SEXP FinalizeCylinderSystem(SEXP ext);
  SEXP CylinderSystem(SEXP R_param, SEXP R_cond);
  SEXP CDigitizeCylinderIntersections(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_delta);
  SEXP IntersectCylinderSystem(SEXP ext, SEXP R_n, SEXP R_z, SEXP R_intern, SEXP R_pl);
  SEXP SimulateCylindersAndIntersect(SEXP R_param, SEXP R_cond, SEXP R_n);

#ifdef __cplusplus
}
#endif


namespace STGM {

class CCylinderSystem {
 public:

  CCylinder::cylinder_type m_type;

  CCylinderSystem(CBox3 &box, double lam, CVector3d &mu, CCylinder::cylinder_type type, int perfect = 1 )
      : m_type(type), m_box(box), m_lam(lam), m_maxR(0), m_mu(mu), num(0), m_perfect(perfect)
  {
    //box.ConstructBoundingPlanes();
  }

  void simSysJoint(R_Calldata d);
  void simCylinderSys(R_Calldata d);
  void simConstCylinderSys(R_Calldata d);
  void simBivariate(R_Calldata d);

  CCylinder simCylinder();

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
