/*
 * IntersectorSpheroidPlane.cpp
 *
 *  Created on: 31.07.2014
 *      Author: franke
 */

#include "Intersector.h"

using namespace std;

namespace STGM
{

  bool Intersector<STGM::CSpheroid>::TestIntersection ()
  {
    int i=0,j=0,k=0;

    //m_spheroid.CalculateMatrixA(m_iF);
    const CMatrix3d &A = m_spheroid.MatrixA();

    /** Get indices according to the given plane for intersections only by XY,XZ,YZ planes  */
    int l;
    for(l=0; l<3; l++)
        if(m_plane.n[l] == 1 || m_plane.n[l] == -1) break;

    switch (l) {
      case 0: // YZ
        i=1;j=2;k=0;
        break;
      case 1: // XZ
        i=0;j=2;k=1;
        break;
      case 2: // XY
        i=0;j=1;k=2;
        break;
    }

    double d = A[i][i]*A[j][j]-SQR(A[i][j]);
    double s[] = { ( A[i][k]*A[j][j]-A[j][k]*A[i][j] ) /d ,
                  ( A[j][k]*A[i][i]-A[i][k]*A[i][j] ) /d };

    double sum = s[0]*s[0]*A[i][i] + s[0]*s[1]*A[i][j] + s[1]*s[0]*A[j][i] + s[1]*s[1]*A[j][j];

    // translate to xj=0 plane
    CVector3d m(m_spheroid.center());
    m[k] = m[k] - m_plane.c;
    double tmp = SQR(m[k])*(A[k][k]-sum);

    //condition for intersection
    return ( tmp <= 1);
 }


  bool Intersector<STGM::CSpheroid>::FindIntersection ()
  {
        int i=0,j=0,k=0;
        const CMatrix3d &A = m_spheroid.MatrixA();

        /** Get indices according to the given plane for intersections only by XY,XZ,YZ planes  */
        int l;
        for(l=0; l<3; l++)
          if(m_plane.n[l] == 1 || m_plane.n[l] == -1) break;

        switch (l) {
          case 0: // YZ
            i=1;j=2;k=0;
            break;
          case 1: // XZ
            i=0;j=2;k=1;
            break;
          case 2: // XY
            i=0;j=1;k=2;
            break;
        }

        double d = A[i][i]*A[j][j]-SQR(A[i][j]);
        double s[] = { ( A[i][k]*A[j][j]-A[j][k]*A[i][j] ) /d ,
                      (  A[j][k]*A[i][i]-A[i][k]*A[i][j] ) /d };

        double sum = s[0]*s[0]*A[i][i] + s[0]*s[1]*A[i][j] + s[1]*s[0]*A[j][i] + s[1]*s[1]*A[j][j];

        // translate to xj=0 plane
        CVector3d m(m_spheroid.center());
        m[k] = m[k] - m_plane.c;
        double tmp = SQR(m[k])*(A[k][k]-sum);

        //condition for intersection
        if( tmp <= 1 )
        {
            CPoint2d m_new;
            CMatrix2d A_new;
            m_new[0]=m[i]+m[k]*s[0];
            m_new[1]=m[j]+m[k]*s[1];
            A_new[0][0]= A[i][i] / (1-tmp);
            A_new[0][1]= A[i][j] / (1-tmp);
            A_new[1][0]= A[j][i] / (1-tmp);
            A_new[1][1]= A[j][j] / (1-tmp);

            /** store ellipse */
            m_ellipse = CEllipse2(A_new,m_new, m_spheroid.Id());

            return true;
       }

      return false;
  }


  bool Intersector<STGM::CSphere>::TestIntersection () {
    double sDist = m_plane.distanceTo(m_sphere.center());
    return fabs(sDist) <= m_sphere.r();
  }

  bool Intersector<STGM::CSphere>::FindIntersection () {
    double sDist = m_plane.distanceTo(m_sphere.center());
    double dist  = fabs(sDist);

    if (dist<=m_sphere.r()) {
      CVector3d center(m_sphere.center());
      center-=sDist*m_plane.n;
      double r = sqrt(fabs(SQR(m_sphere.r())-SQR(dist)));
      m_circle = CCircle3(center,r,m_plane.n,m_sphere.Id());
      return true;
    }
    return false;
  }


  bool Intersector<STGM::CCylinder>::TestIntersection(const CPlane &plane) {
     double UdotN = m_cylinder.u().dot(plane.n);
     double sDist = plane.distanceTo(m_cylinder.origin0());

     if(fabs(UdotN) > 1e-7) // non parallel
     {
           // intersection point
           STGM::CVector3d Ia;
           double t = -sDist/UdotN;
           Ia[0] = m_cylinder.origin0()[0] + t * m_cylinder.u()[0];
           Ia[1] = m_cylinder.origin0()[1] + t * m_cylinder.u()[1];
           Ia[2] = m_cylinder.origin0()[2] + t * m_cylinder.u()[2];

           STGM::CVector3d IaC(Ia);
           IaC -= m_cylinder.center();

           double a = IaC.Length()-0.5*m_cylinder.h();

           if( a <= 0 ) {
               m_type = NON_EMPTY;
               return true;
           }

           // cross product: U x ( N x U)
           STGM::CVector3d w(cross(m_cylinder.u(),STGM::CVector3d(cross(plane.n,m_cylinder.u()))));
           w.Normalize();

           if ( (SQR(a/plane.n.dot(w))-SQR(a)) < SQR(m_cylinder.r())) {
              m_type = NON_EMPTY;
              return true;
           } else {
                 if( fabs(plane.distanceTo(m_cylinder.origin1())) <= m_cylinder.r() || fabs(sDist) <= m_cylinder.r() ) {
                    m_type = NON_EMPTY;
                    return true;
                 }
                 return false;
           }

       } else {
         // parallel
         if( fabs(sDist) <= m_cylinder.r()) {
            m_type = NON_EMPTY;
            return true;
         }
         return false;
       }
   }


  bool Intersector<STGM::CCylinder>::FindIntersection()
   {
     double cosTheta = m_cylinder.u().dot(m_plane.n);
     double absCosTheta = fabs(cosTheta);


     if (absCosTheta > 0.0)
     {
         // The cylinder axis intersects the plane in a unique point.
         if (absCosTheta < 1.0)
         {
             m_type=FindIntersectionType();
             return true;

         } else {
             int id = 1;
             m_type = CIRCLE;
             m_circle1 = STGM::CCircle3(m_cylinder.center(),m_cylinder.r(),m_plane.n,id );
             return true;
         }

      // cylinder parallel to plane
      } else {
          //-> nothing ???
      }
      return false;
   }

  double Intersector<STGM::CCylinder>::GetEllipseSegment(STGM::CVector3d m_center, const STGM::CVector3d &ipt)
  {
    m_center -= ipt;
    double r = sqrt( SQR(m_cylinder.r()) - SQR(m_center.Length()) );

    /** @todo: Get right index here for coordinates according to intersection plane */
    double x = ipt[m_i]+r*m_ellipse.m_minorAxis[m_i];
    double y = ipt[m_j]+r*m_ellipse.m_minorAxis[m_j];

    return acos( (x-m_ellipse.m_center[m_i] + sin(m_cylinder.phi())/cos(m_cylinder.phi())*(y-m_ellipse.m_center[m_j]))
                          /(m_ellipse.m_a*cos(m_cylinder.phi())+m_ellipse.m_a*SQR(sin(m_cylinder.phi()))/cos(m_cylinder.phi())));
  }

  CCircle3 Intersector<STGM::CCylinder>::GetCircle(STGM::CVector3d &spherecenter, double sDist) {
    STGM::CVector3d ctr(spherecenter[0] - sDist*m_plane.n[0],
                        spherecenter[1] - sDist*m_plane.n[1],
                        spherecenter[2] - sDist*m_plane.n[2]);
    double radius = sqrt( SQR(m_cylinder.r()) - SQR(sDist) );
    return STGM::CCircle3(ctr,radius, m_plane.n,1);
  }

  IntersectionType Intersector<STGM::CCylinder>::FindIntersectionType()
  {
     double sDist = m_plane.distanceTo(m_cylinder.origin0());
     double sDist2 = m_plane.distanceTo(m_cylinder.origin1());
     double cosTheta = m_cylinder.u().dot(m_plane.n);
     double absCosTheta = fabs(cosTheta);

     STGM::CVector3d Ia;
     /** intersection point axis with plane: I_a */
     Ia[0] = m_cylinder.origin0()[0] - (sDist/cosTheta)*m_cylinder.u()[0];
     Ia[1] = m_cylinder.origin0()[1] - (sDist/cosTheta)*m_cylinder.u()[1];
     Ia[2] = m_cylinder.origin0()[2] - (sDist/cosTheta)*m_cylinder.u()[2];

     STGM::CVector3d IaC(Ia);
     IaC -= m_cylinder.center();
     double a = IaC.Length()-0.5*m_cylinder.h();

     STGM::CVector3d w(cross(m_cylinder.u(),STGM::CVector3d(cross(m_plane.n,m_cylinder.u()))));
     w.Normalize();
     double NdotW = m_plane.n.dot(w); // cos(theta)

     double b1 = SQR(a/NdotW)-SQR(a);

     m_type = EMPTY;
     m_side = ( fabs(sDist) < fabs(sDist2) ? 1 : -1);

     if(a < 0 || b1 < SQR(m_cylinder.r()))
     {
         m_type=ELLIPSE;
         m_ellipse.m_type = ELLIPSE;
         m_ellipse.m_side = m_side;

         m_ellipse.m_center = Ia;
         m_ellipse.m_majorAxis[0] = m_cylinder.u()[0] - cosTheta*m_plane.n[0];
         m_ellipse.m_majorAxis[1] = m_cylinder.u()[1] - cosTheta*m_plane.n[1];
         m_ellipse.m_majorAxis[2] = m_cylinder.u()[2] - cosTheta*m_plane.n[2];

         m_ellipse.m_minorAxis = cross(m_plane.n,m_ellipse.m_majorAxis); // cross product
         m_ellipse.m_n = m_plane.n;
         m_ellipse.m_a = m_cylinder.r()/absCosTheta;
         m_ellipse.m_b = m_cylinder.r();
         m_ellipse.m_i = m_i;
         m_ellipse.m_j = m_j;

         /** manipulate phi accoriding to the intersecting plane */
         double sroot = sqrt(SQR(m_ellipse.m_majorAxis[m_i])+SQR(m_ellipse.m_majorAxis[m_j]));

         /*///->  TODO: check! why change phi of m_cylinder ? */
         m_cylinder.phi() = ( m_ellipse.m_majorAxis[m_j]>0 ?  acos(m_ellipse.m_majorAxis[m_i]/sroot) : 2*M_PI-acos(m_ellipse.m_majorAxis[m_i]/sroot) );
         m_ellipse.phi() = m_cylinder.phi();
         m_ellipse.m_majorAxis.Normalize();
         m_ellipse.m_minorAxis.Normalize();

         double b2 = SQR((a+m_cylinder.h())/NdotW)-SQR(a+m_cylinder.h());

         if ( b1 < SQR(m_cylinder.r()) && b2 < SQR(m_cylinder.r()) )  {
             /** both ends have circles */
             ipt0[0] = m_ellipse.m_center[0] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[0];
             ipt0[1] = m_ellipse.m_center[1] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[1];
             ipt0[2] = m_ellipse.m_center[2] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[2];

             ipt1[0] = m_ellipse.m_center[0] + m_side*((a+m_cylinder.h())/NdotW)*m_ellipse.m_majorAxis[0];
             ipt1[1] = m_ellipse.m_center[1] + m_side*((a+m_cylinder.h())/NdotW)*m_ellipse.m_majorAxis[1];
             ipt1[2] = m_ellipse.m_center[2] + m_side*((a+m_cylinder.h())/NdotW)*m_ellipse.m_majorAxis[2];

             if(m_side > 0) {
                 m_ellipse.m_psi[0] = GetEllipseSegment(m_cylinder.origin0(),ipt0);
                 m_ellipse.m_psi[1] = GetEllipseSegment(m_cylinder.origin1(),ipt1);
                 m_circle1=GetCircle(m_cylinder.origin0(), sDist);
                 m_circle2=GetCircle(m_cylinder.origin1(), sDist2);
             } else {
                 m_ellipse.m_psi[0] = GetEllipseSegment(m_cylinder.origin1(),ipt0);
                 m_ellipse.m_psi[1] = GetEllipseSegment(m_cylinder.origin0(),ipt1);
                 m_circle2=GetCircle(m_cylinder.origin0(), sDist);
                 m_circle1=GetCircle(m_cylinder.origin1(), sDist2);
             }

             m_type=ELLIPSE_SEGMENT;
             m_ellipse.m_type = ELLIPSE_SEGMENT;
             m_ellipse.m_circle1 = m_circle1;
             m_ellipse.m_circle2 = m_circle2;
             m_ellipse.setReferenceSide();

          } else if(b1 < SQR(m_cylinder.r()) || b2 < SQR(m_cylinder.r())) {
               /** just one end is a circle */
               ipt0[0] = m_ellipse.m_center[0] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[0];
               ipt0[1] = m_ellipse.m_center[1] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[1];
               ipt0[2] = m_ellipse.m_center[2] + m_side*(a/NdotW)*m_ellipse.m_majorAxis[2];

               if(m_side > 0) {
                   m_circle1=GetCircle(m_cylinder.origin0(), sDist);
                   m_ellipse.m_psi[0] = GetEllipseSegment(m_cylinder.origin0(),ipt0);
               }
               else {
                   m_circle1=GetCircle(m_cylinder.origin1(), sDist2);
                   m_ellipse.m_psi[0] = GetEllipseSegment(m_cylinder.origin1(),ipt0);
               }

               m_type=ELLIPSE_ARC;
               m_ellipse.m_type = ELLIPSE_ARC;
               m_ellipse.m_circle1 = m_circle1;
               m_ellipse.setReferenceSide();
         }


     } else {

         /** Ia is lies outside of the cylinder.
          * No intersection with the cylinder
          * but the plane could intersect the caps.
          */

          if (fabs(sDist) <= m_cylinder.r()) {
              // origin0 is the sphere m_center
              m_circle1=GetCircle(m_cylinder.origin0(), sDist);
              m_type=CIRCLE_CAPS;
          }
          if (fabs(sDist2) <= m_cylinder.r()) {
              // origin1 is the sphere m_center
              m_circle1=GetCircle(m_cylinder.origin1(), sDist2);
              m_type=CIRCLE_CAPS;
          }

     }

     return m_type;

  }


} /* namespace STGM */
