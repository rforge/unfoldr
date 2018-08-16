/**
 *  @file  GeometricPrimitives.cpp
 *  @date  02-14-2014
 *
 *  @author: M. Baaske
 */
#define ZERO_TOL 1E-6

#define NDIM 3
#define DISTANCE_VEC(P1,P2,R)                      \
do {                                               \
  int _i;                                          \
  for (_i = 0; _i < NDIM; _i++)                    \
    (R)[_i] = (P1)[_i]-(P2)[_i];                   \
} while(0)

#define Vec3Op(A,assign_op,B,op,C)                 \
(   A[0] assign_op B[0] op C[0],                   \
    A[1] assign_op B[1] op C[1]                    \
)

#define VecDot(A,B)  ((A[0]*B[0]) + (A[1]*B[1]) + (A[2]*B[2]))
#define VecNorm(V)    sqrt(VecDot(V,V))

#define VecScalar(V,assign_op,k)                   \
(   V[0] assign_op k,                              \
    V[1] assign_op k,                              \
    V[2] assign_op k                               \
)

#define MATMULT_MV(V,M,U)                                               \
do {                                                                    \
    int _i, _j;                                                         \
    for (_i = 0; _i < NDIM; _i++) {                                     \
        (V)[_i] = 0.0;                                                  \
        for (_j = 0; _j < NDIM; _j++)                                   \
            (V)[_i] += (M)[_i][_j] * (U)[_j];                           \
    }                                                                   \
} while(0)

#define COPY_R_MATRIX(M,R)                         \
do {                                               \
    int _i, _j;                                    \
    for (_i = 0; _i < NDIM; _i++)                  \
      for (_j = 0; _j < NDIM; _j++)                \
        (M)[_i][_j] = (R)[_j+NDIM*_i];             \
} while(0)

#define sgn(x) ((x) > 0 ? 1 : -1)

#include <R_ext/Lapack.h>
#include "GeometricPrimitives.h"

using namespace std;

static inline double sign(double a,double b) { return a = fabs(a),(b<0)?-a:a; }

/** minimum distance cylinder rods */
void sdm(const double *r12,  const double *u1, const double *u2, const  double *lh1p, const double *lh2p, double *d) {
  double  xmu=0, xla=0,
          lh1=*lh1p, lh2=*lh2p,
          rr = r12[0]*r12[0]+r12[1]*r12[1]+r12[2]*r12[2],
          ru1 = r12[0]*u1[0]+r12[1]*u1[1]+r12[2]*u1[2],
          ru2 = r12[0]*u2[0]+r12[1]*u2[1]+r12[2]*u2[2],
          u1u2 = u1[0]*u2[0]+u1[1]*u2[1]+u1[2]*u2[2],
          cc = 1.0-SQR(u1u2);

  // Checking whether the rods are or not parallel:
  // The original code is modified to have symmetry:

   if(cc<ZERO_TOL) {
    if(ru1!=0 && ru2!=0) {
        xla= ru1/2;
        xmu= -ru2/2;
    } else {
        xla=0;
        xmu=0;
    }
   } else {
    // Step 1
    xla= (ru1-u1u2*ru2)/cc;
    xmu= (-ru2+u1u2*ru1)/cc;
   }

  // Step 2
  if( fabs(xla)>lh1 || fabs(xmu)>lh2 ) {
  // Step 3 - 7
    if(fabs(xla)-lh1>fabs(xmu)-lh2) {
     xla= sign(lh1,xla);
     xmu= xla*u1u2-ru2;
     if( fabs(xmu)>lh2 ) xmu= sign(lh2,xmu);
    }
    else {
     xmu= sign(lh2,xmu);
     xla= xmu*u1u2+ru1;
     if( fabs(xla)>lh1 ) xla= sign(lh1,xla);
    }
   }
   // Step 8
   *d=rr+SQR(xla)+SQR(xmu)-2*xla*xmu*u1u2 + 2*xmu*ru2 - 2*xla*ru1;
}


double solveQ(double p, double q) {
  double r=0, D=SQR(0.5*p)-q;
  if(D>0) {
     r = -0.5*p-sgn(p)*sqrt(D);
     return ( r<q/r ? q/r: r);
  }
  return -1;
}

double contactRadius(double *t, double *rr, double li, double lj, double ri, double rj) {
  double z[3] = {0.0,0.0,1.0}, a = ri+rj,
          tmp = 0, x = 0, b = 0, ti = 0, tj = 0,
         rmax = 0;                                                       /* final contact radius*/

  double r = VecNorm(rr);
  VecScalar(rr,/=,r);

  double zt = VecDot(z,t),
         rz = VecDot(rr,z),
         rt = VecDot(rr,t);

  if(zt<1) {
      double Ai = (rz-rt*zt)/(1-SQR(zt));
      double Aj = (rz*zt-rt)/(1-SQR(zt));
      tmp = a/sqrt(1+SQR(Ai)+SQR(Aj)+2*(Aj*rt-Ai*rz)-2*Ai*Aj*zt);
      ti = tmp*Ai;
      tj = tmp*Aj;

      if( (ti>-li && ti<li) && (tj>-lj && tj < lj) )
        rmax = tmp;
  }

  if(rz<1) {
      x = 2*lj*(rt-rz*zt)/(1-SQR(rz));
      b = (SQR(lj)*SQR(zt)+SQR(lj)-2*SQR(lj)*SQR(zt)-SQR(a))/(1-SQR(rz));
      tmp = solveQ(x,b);
      if(tmp>0) {
          ti = tmp*rz+lj*zt;
          if( ti>-li && ti<li && tmp>rmax)
            rmax = tmp;
      }

      ///x = 2*lj*(rz*zt-rt)/(1-SQR(rz));
      tmp = solveQ(-x,b);
      if(tmp>0) {
          ti = tmp*rz-lj*zt;
          if( ti>-li && ti<li && tmp>rmax)
            rmax = tmp;
      }

  }

  if(rt<1) {
      x = 2*li*(rt*zt-rz)/(1-SQR(rt));
      b = (SQR(li)*SQR(zt)+SQR(li)-2*SQR(li)*SQR(zt)-SQR(a))/(1-SQR(rt));
      tmp = solveQ(x,b);
      if(tmp>0) {
          tj = li*zt-tmp*rt;
          if( tj>-lj && tj<lj && tmp>rmax)
            rmax = tmp;
      }

      ///x = 2*li*(rz-rt*zt)/(1-SQR(rt));
      tmp = solveQ(-x,b);
      if(tmp>0) {
         tj = -li*zt-tmp*rt;
         if( tj>-lj && tj<lj && tmp>rmax)
           rmax = tmp;
      }

  }

  /* do not change order! */
  x = 2*(lj*rt-li*rz);
  b = SQR(li)+SQR(lj)-2*li*lj*zt-SQR(a);
  tmp = solveQ(x,b);
  if(tmp>0 && tmp>rmax)
    rmax = tmp;

  ///x = 2*(-lj*rt+li*rz);
  tmp = solveQ(-x,b);
  if(tmp>0 && tmp>rmax)
    rmax = tmp;

  x = 2*(-lj*rt-li*rz);
  b = SQR(li)+SQR(lj)+2*li*lj*zt-SQR(a);
  tmp = solveQ(x,b);
  if(tmp>0 && tmp>rmax)
    rmax = tmp;

  ///x = 2*(lj*rt+li*rz);
  tmp = solveQ(-x,b);
  if(tmp>0 && tmp>rmax)
    rmax = tmp;

  return rmax;
}

/** cylinder contact radius */
void ContactRadius(double *u, double *li, double *lj, double *ri, double *rj,double *R, double *d, double *rmax)
{
  double m[NDIM][NDIM];
  COPY_R_MATRIX(m,R);

  double rr[3],t[3];
  MATMULT_MV(t,m,u);
  MATMULT_MV(rr,m,d);

  *rmax=contactRadius(t,rr, *li, *lj, *ri, *rj);
}


namespace STGM {

  typedef const double *array_t;
  typedef const double (*matrix_t)[3];

  void real_eval(double *a, int *n, double *evalf, int *err) {
	  // size(evalf): n
	  // result: evec = a

	  int lda = *n,  lwork = 3*lda-1;
	  double *work = Calloc(lwork, double);
	  F77_NAME(dsyev)("V","U", &lda, a, &lda, evalf, work, &lwork, err);
	  Free(work);

  }

  double CWindow::PointInWindow(STGM::CVector2d point)  {
     STGM::CVector2d diff(point);
     diff -= m_center;

     double sqrDistance = 0, delta=0;
     STGM::CVector2d closest;

     for (int i = 0; i < 2; ++i)
     {
         closest[i] = diff.dot( *(m_axis[i]));
         if (closest[i] < -m_extent[i])
         {
             delta = closest[i] + m_extent[i];
             sqrDistance += delta*delta;
             closest[i] = -m_extent[i];
         }
         else if (closest[i] > m_extent[i])
         {
             delta = closest[i] - m_extent[i];
             sqrDistance += delta*delta;
             closest[i] = m_extent[i];
         }
     }

     STGM::CVector2d mClosestPoint = m_center;
     for (int i = 0; i < 2; ++i)
       mClosestPoint += closest[i]* *(m_axis[i]);

     return sqrDistance;
  }


  CVector3d PerpendicularVector(const CVector3d &a) {
      if (fabs(a[2]) > fabs(a[0]))
        return CVector3d(0.0, -a[2], a[1]);
      else
        return CVector3d(-a[1], a[0], 0.0);
  }

  CMatrix3d RotationMatrixFrom001(CVector3d v) {
    v.Normalize();

    CVector3d o1(PerpendicularVector(v));
    o1.Normalize();

    CVector3d o2(cross(v,o1));
    o2.Normalize();

    CMatrix3d r;

    r[0][0] = o1[0];
    r[0][1] = o1[1];
    r[0][2] = o1[2];

    r[1][0] = o2[0];
    r[1][1] = o2[1];
    r[1][2] = o2[2];

    r[2][0] = v[0];
    r[2][1] = v[1];
    r[2][2] = v[2];

    return r;
  }

  double CSpheroid::spheroidDistanceAsCylinder(CSpheroid &sp) const {
      STGM::CVector3d w(0,0,0);
      double d=0, lh1=1.0e-7,lh2=1.0e-7;

      DISTANCE_VEC(m_center.ptr(),sp.center().ptr(),w.ptr());

      if(m_crack && sp.m_crack)  {  // E_i, E_j
         lh1 = m_b-m_a;
         lh2 = sp.b()-sp.a();
         sdm(w.ptr(),m_u.ptr(),sp.u().ptr(),&lh1,&lh2,&d);
         return MAX(sqrt(d) - m_a - sp.a(),0.0);

      } else if(!m_crack && sp.m_crack)  {  // D_i, E_j
         lh2 = sp.b()-sp.a();
         sdm(w.ptr(),m_u.ptr(),sp.u().ptr(),&lh1,&lh2,&d);
         return MAX(sqrt(d) - sp.a(),0.0);

      } else if(m_crack && !sp.m_crack)  {  // E_i, D_j
          lh1 = m_b-m_a;
          sdm(w.ptr(),m_u.ptr(),sp.u().ptr(),&lh1,&lh2,&d);
          return MAX(sqrt(d) - m_a,0.0);

      } else  {//if(!m_crack && !sp.m_crack)  // D_i, D_j
          sdm(w.ptr(),m_u.ptr(),sp.u().ptr(),&lh1,&lh2,&d);
          return MAX(sqrt(d),0.0);
      }

    }


  void CSpheroid::ComputeMatrixA() {
     m_A.nullify();
	 //m_A[0][0] = m_A[1][1] = 1.0 / SQR(m_a);
     m_A[0][0] = 1.0 / SQR(m_a);
     m_A[1][1] = 1.0 / SQR(m_c);
     m_A[2][2] = 1.0 / SQR(m_b);

     CMatrix3d R = RotationMatrixFrom001(m_u);

     m_A = m_A * R;
     R.Transpose();
     m_A = R * m_A;
   }

   PointVector2d CCircle3::getMinMaxPoints()  {
     PointVector2d p;
     p.push_back(getMinMax_X());
     p.push_back(getMinMax_Y());
     return p;
   }

   PointVector2d CCircle3::getExtremePointsCircle() {
     PointVector2d p;
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[0]));
     p.push_back(CPoint2d(m_center[m_i],getMinMax_Y()[1]));
     p.push_back(CPoint2d(getMinMax_X()[0],m_center[m_j]));
     p.push_back(CPoint2d(getMinMax_X()[1],m_center[m_j]));
     return p;
   }


   double CCylinder::cylinderDistance(CCylinder &sp) const {
       double d=0, lh1=.5*m_h, lh2=.5*sp.m_h;

       STGM::CVector3d w(0,0,0), z(0,0,1);
       DISTANCE_VEC(m_center.ptr(),sp.center().ptr(),w.ptr());

       sdm(w.ptr(),m_u.ptr(),sp.u().ptr(),&lh1,&lh2,&d);
       return MAX(sqrt(d)-m_r-sp.m_r,0.0);
   }

   CEllipse2 CCylinder::crackProjection() const {
	   STGM::CVector2d ctr(m_center[0],m_center[1]);
       STGM::CVector3d z(cos(m_phi)*m_u[2], sin(m_phi)*m_u[2], sin(-m_phi)*m_u[1]-cos(-m_phi)*m_u[0]);
       STGM::CVector3d pz(z), pv(cross(m_u,z));    // minor, major

       pz.Normalize();
       pv.Normalize();

       pz *= m_r;
       pz += m_center;
       pv *= m_r;
       pv += m_center;

       STGM::CVector2d axis1(pz[0]-ctr[0],pz[1]-ctr[1]);
       STGM::CVector2d axis2(pv[0]-ctr[0],pv[1]-ctr[1]);

       double zlen = axis1.Length();
       double vlen = axis2.Length();

       axis1.Normalize();
       axis2.Normalize();

       return STGM::CEllipse2(ctr,axis2,axis1,zlen,vlen,m_id);
   }


   double CCylinder::delamProjection(PointVector2d &P, int npoints) {
      int np = std::floor(MAX(8.0,(double) npoints-4.0)/2.0);
      CVector3d n(0,0,1), u(cross(m_u,n));

      /// projected rectangle points
      u.Normalize();
      P.push_back(CPoint2d(m_center[0]+m_r*u[0]+0.5*m_h*m_u[0],m_center[1]+m_r*u[1]+0.5*m_h*m_u[1]));
      P.push_back(CPoint2d(m_center[0]-m_r*u[0]+0.5*m_h*m_u[0],m_center[1]-m_r*u[1]+0.5*m_h*m_u[1]));
      P.push_back(CPoint2d(m_center[0]+m_r*u[0]-0.5*m_h*m_u[0],m_center[1]+m_r*u[1]-0.5*m_h*m_u[1]));
      P.push_back(CPoint2d(m_center[0]-m_r*u[0]-0.5*m_h*m_u[0],m_center[1]-m_r*u[1]-0.5*m_h*m_u[1]));

      /// rectangle sides
      CPoint2d A,B;
      Vec3Op(A,=,P[0],-,P[2]);
      Vec3Op(B,=,P[0],-,P[1]);

      CCircle3 c1(m_origin0,m_r,n);
      c1.samplePoints(P,np);
      CCircle3 c2(m_origin1,m_r,n);
      c2.samplePoints(P,np);

      /// projection area: 2*area of the half circles times
      /// the area of the rectangle in between the closing caps
      return c1.area() + A.Length()*B.Length();

   }

   double CCylinder::projectedPointsWithArea(PointVector2d &P, int npoints) {
     /// delam, full projection
     if(m_crack) {
       /// for spheres as spherocylinders
       if(!std::strcmp(m_label,"F")) {
           CCircle3 circle(m_center,m_r);
           circle.samplePoints(P,npoints);
           return M_PI*SQR(m_r);
       }
       /// only for full spherocylinders
       return delamProjection(P,npoints);
     /// crack, circle projection (same as for spheroids)
     } else {
       CEllipse2 e = crackProjection();
       e.samplePoints(P,npoints);
       return e.area();
     }
   }


  void CBox3::setExtent(double a, double b, double c) {
     m_extent[0] = 0.5*a;
     m_extent[1] = 0.5*b;
     m_extent[2] = 0.5*c;

     m_size[0] = a;
     m_size[1] = b;
     m_size[2] = c;
  }

  void CBox3::ConstructBoundingPlanes()
  {
    /** sorted order of normal vectors (planes) -  do not change! */
	  m_planes.push_back(CPlane(STGM::CVector3d(1,0,0),m_low[0])); // left
	  m_planes.push_back(CPlane(STGM::CVector3d(1,0,0),m_up[0]));  // right in x direction
	  m_planes.push_back(CPlane(STGM::CVector3d(0,1,0),m_low[1])); // front
	  m_planes.push_back(CPlane(STGM::CVector3d(0,1,0),m_up[1]));  // back in y direction
	  m_planes.push_back(CPlane(STGM::CVector3d(0,0,1),m_low[2])); // bottom
	  m_planes.push_back(CPlane(STGM::CVector3d(0,0,1),m_up[2]));  // top
  }

  void CBox3::ConstructBoxLateralPlanes()
  {
	  /** sorted order of normal vectors (planes) -  do not change! */
	  m_lateral_planes.push_back(CPlane(STGM::CVector3d(1,0,0),m_low[0])); // left
	  m_lateral_planes.push_back(CPlane(STGM::CVector3d(1,0,0),m_up[0]));  // right in x direction
	  m_lateral_planes.push_back(CPlane(STGM::CVector3d(0,1,0),m_low[1])); // front
	  m_lateral_planes.push_back(CPlane(STGM::CVector3d(0,1,0),m_up[1]));  // back in y direction
  }

  STGM::CEllipse2 CSpheroid::spheroidProjection() const {
    if(m_crack)
      return delamProjection();
    else
      return crackProjection();
  }


  STGM::CEllipse2 CSpheroid::delamProjection() const
  {
    STGM::CMatrix2d A_new;
    A_new[0][0]= m_A[0][0];
    A_new[0][1]= m_A[0][1];
    A_new[1][0]= m_A[1][0];
    A_new[1][1]= m_A[1][1];

    array_t m = m_center.ptr();
    STGM::CVector2d center(m[0],m[1]);

    // TODO: check! phi.
    //return STGM::CEllipse2(A_new, center, m_id, 0.5*M_PI);
    return STGM::CEllipse2(A_new, center, m_id);
  }

  STGM::CEllipse2 crackProjection(STGM::CVector3d &center, STGM::CVector3d &u, double a, double phi, int id)
  {
	STGM::CVector2d ctr(center[0],center[1]);
	STGM::CVector3d z(cos(phi)*u[2], sin(phi)*u[2], sin(-phi)*u[1]-cos(-phi)*u[0]);
    STGM::CVector3d pz(z), pv(cross(u,z));  // minor and major
    pz.Normalize();
    pv.Normalize();
    pz *= a;
    pz += center;
    pv *= a;
    pv += center;

    STGM::CVector2d axis1(pz[0]-ctr[0],pz[1]-ctr[1]);
    STGM::CVector2d axis2(pv[0]-ctr[0],pv[1]-ctr[1]);

    double zlen = axis1.Length();
    double vlen = axis2.Length();
    axis1.Normalize();
    axis2.Normalize();

    return STGM::CEllipse2(ctr,axis2,axis1,zlen,vlen,id);
  }

  STGM::CEllipse2 CSpheroid::crackProjection() const
  {
	    STGM::CVector2d ctr(m_center[0],m_center[1]);
	    STGM::CVector3d z(cos(m_phi)*m_u[2],
                          sin(m_phi)*m_u[2],
                          sin(-m_phi)*m_u[1]-cos(-m_phi)*m_u[0]);

        STGM::CVector3d pz(z),               // minor
                        pv(cross(m_u,z));    // major

        pz.Normalize();
        pv.Normalize();

        pz *= m_a;
        pz += m_center;
        pv *= m_a;
        pv += m_center;

        STGM::CVector2d axis1(pz[0]-ctr[0],pz[1]-ctr[1]);
        STGM::CVector2d axis2(pv[0]-ctr[0],pv[1]-ctr[1]);

        double zlen = axis1.Length();
        double vlen = axis2.Length();

        axis1.Normalize();
        axis2.Normalize();

        return STGM::CEllipse2(ctr,axis2,axis1,zlen,vlen,m_id);
  }

  bool CEllipse3::isInside(double x, double y) {
       if(m_type == 7) {
          return isInsideEllipse(x,y);
       } else if(m_type == 8) {
            int side1 = whichSide(STGM::CPoint2d(x,y),0);
            if( side1 == m_side0 || side1 == 0)
              return isInsideEllipse(x,y);
            else
              return m_circle1.isInside(x,y);
       } else if(m_type==9) {
           int side1 = whichSide(STGM::CPoint2d(x,y),0);
           int side2 = whichSide(STGM::CPoint2d(x,y),1);
           if( (side1 == m_side0 || side1 == 0) && (-side2 == m_side0 || side2 == 0) ) {
             return isInsideEllipse(x,y);
           } else {
             return ( m_circle1.isInside(x,y) || m_circle2.isInside(x,y));
           }

       }
       return false;
  }


  bool CEllipse3::isInWindow(STGM::CWindow &win) {
      std::vector<STGM::CPoint2d> p = getMinMaxPoints();

      if( (win.PointInWindow( STGM::CVector2d(p[0][0],p[1][0]) ) == 0) &&
          (win.PointInWindow( STGM::CVector2d(p[0][0],p[1][1]) ) == 0) &&
          (win.PointInWindow( STGM::CVector2d(p[0][1],p[1][0]) ) == 0) &&
          (win.PointInWindow( STGM::CVector2d(p[0][1],p[1][1]) ) == 0))
      {
        return true;
      }
      return false;
    }

    PointVector2d CEllipse3::getEllipseExtremePoints()
    {
         //Rprintf("getEllipseExtremePoints...\n");
         std::vector<STGM::CPoint2d> p;
         p.push_back(getMaxEllipsePoint_X());
         p.push_back(getMinEllipsePoint_X());
         p.push_back(getMaxEllipsePoint_Y());
         p.push_back(getMinEllipsePoint_Y());

  #if DEBUG
         for(int i=0; i<p.size();i++)
           Rprintf("[ %f %f ], \n",p[i][0], p[i][1]);
         Rprintf("\n");

  #endif
         return p;
     }


    /** @brief Minimum and maximum coordinates
     *
     * @return
     */
     std::vector<STGM::CPoint2d> CEllipse3::getMinMaxPoints()
     {
         //Rprintf("getMinMaxPoints... \n");
         std::vector<STGM::CPoint2d> p;
         p.reserve(2);

         /** ELLIPSE */
         if(m_type == 7) {
             p.push_back(getMinMax_X());
             p.push_back(getMinMax_Y());

  #if DEBUG
             Rprintf("return p...\n");
             for(int i=0; i<p.size();i++)
                Rprintf("[ %f %f ], \n",p[i][0], p[i][1]);
             Rprintf("\n");
  #endif
             return p;

         /** ELLIPSE_ARC and ELLIPSE_SEG*/
         } else if(m_type==8 || m_type==9) {
             int side1=-2, side2=-2;
             /** Ellipse points*/
             std::vector<STGM::CPoint2d> py; /** at end y sorted coordinates */
             std::vector<STGM::CPoint2d> pp = getEllipseExtremePoints();

             for(PointIterator it = pp.begin(); it != pp.end(); ++it) {
                 side1 = whichSide(*it,0);
                 if(m_type==8) {
                     /** point must not lie in the circle cap*/
                     if( side1 == m_side0 || side1 == 0)
                       py.push_back(*it);
                 } else if(m_type==9) {
                     side2 = whichSide(*it,1);
                     /** point must not lie in both circle caps */
                     if( (side1 == m_side0 || side1 == 0) && ( -side2 == m_side0 || side2 == 0) )
                         py.push_back(*it);
                 }
                 //Rprintf("side0 %d == side1 %d, side2 %d\n", m_side0,side1,side2);
             }

  #if DEBUG
             Rprintf("\n");
             for(int i=0; i<py.size();i++)
                Rprintf("[ %f %f ], \n",py[i][0], py[i][1]);
              Rprintf("\n");
  #endif
             /** Circle points*/
             std::vector<STGM::CPoint2d> q = m_circle1.getExtremePointsCircle();

  #if DEBUG
             Rprintf("\n");
             for(int i=0; i<q.size();i++)
                Rprintf("[ %f %f ], \n",q[i][0], q[i][1]);
              Rprintf("\n");
  #endif

             for(PointIterator it = q.begin(); it != q.end(); ++it)
               py.push_back(*it);

             /** push extreme points of the second circle cap */
             if(m_type==9) {
                q = m_circle2.getExtremePointsCircle();
                for(PointIterator it = q.begin(); it != q.end(); ++it)
                  py.push_back(*it);
             }

             PointVector2d::iterator max_x = std::max_element(py.begin(), py.end(), compareX);
             PointVector2d::iterator min_x = std::min_element(py.begin(), py.end(), compareX);
             PointVector2d::iterator max_y = std::max_element(py.begin(), py.end(), compareY);
             PointVector2d::iterator min_y = std::min_element(py.begin(), py.end(), compareY);

  #if DEBUG
             Rprintf("sort y...\n");
             for(int i=0; i<py.size();i++)
               Rprintf("[ %f %f ], \n",py[i][0], py[i][1]);
             Rprintf("\n");

  #endif
             p.push_back( STGM::CPoint2d( (*min_x)[0],(*max_x)[0]) ); // x coordinate sorted as [min, max]
             p.push_back( STGM::CPoint2d( (*min_y)[1],(*max_y)[1]) ); // y coordinate sorted as [min, max]

         }
         return p;
    }


} /* STGM */
