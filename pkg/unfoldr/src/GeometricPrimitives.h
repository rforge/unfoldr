/**
 * @file GeomoetricPrimitives.h
 * @date 2014-02-11
 *
 * @brief classes for geometric objects with methods for
 *        intersection, projection and certain distance definitions
 *
 * @author: M. Baaske
 */

#ifndef GEOMETRIC_PRIMITIVES_H_
#define GEOMETRIC_PRIMITIVES_H_

#include "Vector.h"

extern "C" void sdm(const double *r12,  const double *u1, const double *u2, const  double *lh1p, const double *lh2p, double *d);
extern "C" void ContactRadius(double *u, double *li, double *lj, double *ri, double *rj,double *R, double *d, double *rmax);

namespace STGM {

// Information about the intersection set
  enum IntersectionType  {
    EMPTY=0,          // 0
    NON_EMPTY,        // 1
    POINT,            // 2
    LINE,             // 3
    LINES,            // 4
    DISC,       	  // 5
    CAP,		      // 6
    ELLIPSE,          // 7
    ELLIPSE_ARC,      // 8
    ELLIPSE_SEGMENT,  // 9
    ELLIPSE_2D        // 10
  };

  /** some type definitions */
  typedef std::vector<STGM::CPoint2d> PointVector2d;
  typedef std::vector<STGM::CPoint2d>::iterator PointIterator;

  void real_eval(double *a, int *n, double *evalf, int *err);

  /**
   * Calulate rotation matrix which
   * transforms (0,0,1) into vector v
   */
  CMatrix3d RotationMatrixFrom001(CVector3d v);

  class CWindow
   {
   public:
     /**
      * @brief Constructor
      *        A window within the first quadrant
      *
      * @param a,b Length of each direction u=(1,0) and v=(0,1)
      */

     CWindow(double a, double b)
      : m_size(STGM::CPoint2d(a,b)),
		m_center(STGM::CVector2d(0.5*a,0.5*b)),
		m_low(0,0),
		m_up(a,b),
		m_u(STGM::CVector2d(1.0,0.0)),
        m_v(STGM::CVector2d(0.0,1.0))
     {
         m_axis[0] = &m_u;
         m_axis[1] = &m_v;
         m_extent[0] = 0.5*a;
         m_extent[1] = 0.5*b;
     }

     CWindow(STGM::CPoint2d p)
       : m_size(p),
		 m_center(STGM::CVector2d(0.5*p[0],0.5*p[1])),
		 m_low(0,0),
		 m_up(p[0],p[1]),
		 m_u(STGM::CVector2d(1.0,0.0)),
		 m_v(STGM::CVector2d(0.0,1.0))
	  {
		  m_axis[0] = &m_u;
		  m_axis[1] = &m_v;
		  m_extent[0] = 0.5*p[0];
		  m_extent[1] = 0.5*p[1];
	  }


     CWindow(double xrange[2], double yrange[2])
      :  m_size(STGM::CPoint2d(std::fabs(xrange[1]-xrange[0]),std::fabs(yrange[1]-yrange[0]))),
		 m_center(STGM::CVector2d(xrange[1]-0.5*m_size[0],yrange[1]-0.5*m_size[1])),
         m_low(xrange[0],yrange[0]),
         m_up(xrange[1],yrange[1]),
		 m_u(STGM::CVector2d(1.0,0.0)),
		 m_v(STGM::CVector2d(0.0,1.0))
     {
    	 m_axis[0] = &m_u;
    	 m_axis[1] = &m_v;
    	 m_extent[0] = 0.5*m_size[0];
    	 m_extent[1] = 0.5*m_size[1];
     }

     /**
      * @brief Constructor
      *
      * @param center The center of the window
      * @param u      Oriented first axis , e.g. x
      * @param v      Oriented second axis, e.g. y
      * @param a,b    Length of each direction u and v
      */

     CWindow(double center[2], double a, double b )
       : m_size(STGM::CPoint2d(a,b)),
		 m_center(STGM::CVector2d(center[0],center[1])),
		 m_low(center[0]-0.5*a,center[1]-0.5*b),
		 m_up(center[0]+0.5*a,center[1]+0.5*b),
		 m_u(STGM::CVector2d(1.0,0.0)),
         m_v(STGM::CVector2d(0.0,1.0))
     {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_extent[0] = 0.5*a;
        m_extent[1] = 0.5*b;
     }

     CWindow(const CWindow &other) :
    	 m_size(other.m_size),
    	 m_center(other.m_center),
		 m_low(other.m_low),
		 m_up(other.m_up),
		 m_u(other.m_u), m_v(other.m_v)
     {
       m_axis[0] = &m_u;
       m_axis[1] = &m_v;
       m_extent[0] = other.m_extent[0];
       m_extent[1] = other.m_extent[1];
     }

     /** Copy Assignment Operator */
     CWindow& operator= (const CWindow& other) {
       if (this != &other) {
         CWindow tmp(other); // re-use copy-constructor
         *this = tmp;
       }
       return *this;
     }

     virtual ~CWindow() {};

     double PointInWindow(STGM::CVector2d point);

     STGM::CPoint2d  m_size;
     STGM::CVector2d m_center;
     STGM::CVector2d m_low, m_up;
     STGM::CVector2d m_u, m_v, *m_axis[2];
     double m_extent[2];

   };

  /**
   * @brief A base class
   *
   */

  class CGeometry
  {
  public:
    CGeometry() {};
    virtual ~CGeometry() {};

    /** not really meaningful, but overloaded anyway! */
    virtual bool isInWindow(CWindow &win)
     { return win.PointInWindow(CVector2d(0,0)) == 0; };

    virtual bool isInside(double x=0.0, double y=0.0)
     { return x != y; };

    virtual void move(CVector2d &v)
     { v[0]=0; v[1]=0;
       return;
     }								/* move 2D object relative to window [0,0] */

    virtual PointVector2d getMinMaxPoints() { return PointVector2d(); }
  };

   /**
    *  @brief A simple Bounding rectangle aligned to
    *         x and y axis
    */

   class CBoundingRectangle : public CGeometry
   {
     public:
       CBoundingRectangle() :
         m_ymin(0), m_ymax(0), m_xmin(0), m_xmax(0)
       {
       }

       CBoundingRectangle(int ymin, int ymax, int xmin, int  xmax) :
           m_ymin(ymin),m_ymax(ymax), m_xmin(xmin), m_xmax(xmax)
       {
       };

       virtual ~CBoundingRectangle() {};

       int m_ymin, m_ymax;
       int m_xmin, m_xmax;
   };


   /**
    *  @brief Plane
    */
   class CPlane : public CGeometry
   {
   public:

       CPlane () :  n( CVector3d(0,0,1) ), c(0) {}
       virtual ~CPlane () {};

       CPlane ( const CVector3d &_n, const CVector3d &_p) : n(_n)
       {
          c=n.dot(_p);
       }

       CPlane (const CVector3d &_n, const double &_c = 0.0) : n(_n), c(_c)
       {}

       inline double distanceTo ( const CVector3d &p) const {
         return n.dot(p)-c;
       }

       /**
        * @param p point on plane
        * @return -1,0,1
        */
        inline int whichSide ( const CVector3d &p) const {
          double dist = distanceTo(p);
            if (dist < 0)
                return -1; // p on negative side
            else if (dist > 0)
                return +1; // p on positive side
            else
                return 0;  // p on the plane
        }

        /**
         *
         * @return index of '1'
         */
        inline int idx() const {
           int l=0;
           for(l=0; l<3; l++) if(n[l] == 1 || n[l] == -1) break;
           return l;
        }


        void getPlaneIdx(int &i, int &j) {
			switch(idx()) {
			  case 0: i=1; j=2; break; // YZ
			  case 1: i=0; j=2; break; // XZ
			  case 2: i=0; j=1; break; // XY
			}
		}

        void getPlaneIdx(int &i, int &j, int &k) {
			switch(idx()) {
			  case 0: i=1; j=2; k=0; break; // YZ
			  case 1: i=0; j=2; k=1; break; // XZ
			  case 2: i=0; j=1; k=2; break; // XY
			}
		}

        CVector3d n;
        double c;
   };

   typedef std::vector<CPlane> LateralPlanes;

   /**
    * @brief Circle 3d
    */
   class CCircle3 : public CGeometry
   {
   public:
     CCircle3 () :
       m_center(CVector3d(0,0,0)), m_n(CVector3d(0,0,1)), m_plane(CPlane()), m_radius(0), m_id(0)
     {
       setPlaneIdx();
     };

     ~CCircle3 (){};

     CCircle3(CVector3d &_center, double _radius, CVector3d &_n, int id = 0)
       : m_center(_center), m_n(_n), m_plane(CPlane(_n)), m_radius(_radius), m_id(id)
     {
       setPlaneIdx();
     }

     CCircle3(CVector3d &_center, double _radius )
       : m_center(_center), m_n(CVector3d(0,0,1)), m_plane(CPlane()), m_radius(_radius), m_id(0)
      {
        setPlaneIdx();
      }

     inline void setPlaneIdx() {
        switch(m_plane.idx()) {
          case 0: m_i=1; m_j=2; break; // YZ
          case 1: m_i=0; m_j=2; break; // XZ
          case 2: m_i=0; m_j=1; break; // XY
        }
     }

     void setPlane(CPlane& plane) {
       m_plane = plane;
       m_n = m_plane.n;
       setPlaneIdx();
     }

     CPoint2d PointOnCircle(const double t) {
       return CPoint2d( m_center[m_i]+m_radius*cos(t), m_center[m_j]+m_radius*sin(t) );
     }

     void move(CVector2d &v) {
    	 m_center[m_i] -= v[0];
    	 m_center[m_j] -= v[1];
     }

     /**
      * @brief   sample the ellipse border with np points
      *
      * @param P         Point vector of sampled points
      * @param np        number of points to sample
      */
     void samplePoints(PointVector2d &P, int np) {
       double t=0.0, s=2.0*M_PI / (double)np;
       for(int k=0; k<np; ++k) {
           P.push_back( PointOnCircle(t) );
           t += s;
       }
     }

     inline int Id() { return m_id; }

     CVector3d & center() { return m_center; }
     const CVector3d center() const { return m_center; }

     double & r() { return m_radius; }
     double r() const { return m_radius; }

     CVector3d& n() { return m_n; }
     const CVector3d n() const { return m_n; }

     inline CPoint2d getMinMax_X() { return CPoint2d(m_center[m_i]-m_radius,m_center[m_i]+m_radius);  }
     inline CPoint2d getMinMax_Y() { return CPoint2d(m_center[m_j]-m_radius,m_center[m_j]+m_radius);  }

     PointVector2d getMinMaxPoints();
     PointVector2d getExtremePointsCircle();

     inline bool isInside(double x, double y)  { return  SQR(x-m_center[m_i])+SQR(y-m_center[m_j])<=SQR(m_radius); };

     CBoundingRectangle & getBoundingRectangle() { return m_br; };

     inline bool isInWindow(CWindow &win) {
      if( win.PointInWindow(CVector2d(m_center[m_i],m_center[m_j])) == 0)
       return true;
      return false;
     }

     inline double area() { return M_PI*m_radius*m_radius; }

     CVector3d m_center, m_n;
     CPlane m_plane;
     double m_radius;
     int m_i, m_j;
     CBoundingRectangle m_br;
     int m_id;


   };


   class CSphere :  public CGeometry {
      public:
        const char *m_label;

        CSphere(double _x, double _y, double _z, double _r, int id = 0, const char *label = "N", int interior = 1)
         : m_label(label), m_crack(0), m_id(id), m_center(CVector3d(_x,_y,_z)),
		   m_r(_r), m_interior(interior)
        {
        }

        CSphere(const CVector3d &_center, const double _r, const int id, const char *label = "N", int interior=1)
         : m_label(label), m_crack(0), m_id(id), m_center(_center),
		   m_r(_r), m_interior(interior)
        {
        }

        ~CSphere() {};

        CVector3d& center() { return m_center; }
        const CVector3d& center() const { return m_center; }

        inline int Id() { return m_id; }

        double& r() { return m_r; }
        const double& r() const { return m_r; }

        int& interior() { return m_interior; }
        const int& interior() const { return m_interior; }

        const char * label() const { return m_label; }
        void setCrackType(int crack) { m_crack=crack; }

        /**
         * @brief "crack" is alway "delam"
         *
         * @param n  normal vector of plane to project on
         * @return   CCircle3
         */
        CCircle3 project(CVector3d &n) { return CCircle3(m_center,m_r,n);  }
        inline double projectionArea() { return M_PI*SQR(m_r); }

        double sphereDistance(CSphere &s) const {
          CVector3d d(m_center);
          d -= s.m_center;
          return d.Length();
        }

        inline double distance(CSphere &s) {  return sphereDistance(s); }

        /**
          * @brief Sample points and return area of projection
          *
          * @param P       Point vector of points
          * @param np      number of points to sample
          * @return        area of ellipse
          */
        inline double projectedPointsWithArea(PointVector2d &P, int np) {
          CVector3d n(0,0,1);
          CCircle3 circle(m_center,m_r,n);
          circle.samplePoints(P,np);
          return M_PI*SQR(m_r);
        }


      private:
        int m_crack, m_id;
        CVector3d m_center;
        double m_r;
        int m_interior;

      };

   typedef std::vector<CSphere> Spheres;

  /**
   * @brief Ellipse in 2d
   */

  class CEllipse2 : public CGeometry
  {
  public:

    CEllipse2() :
      m_center(STGM::CVector2d()), m_a(0), m_b(0),  m_phi(0), m_id(0), m_type(10)
    {
    };

    virtual ~CEllipse2() {};

    /**
     * @brief  Constructor, calculate major and minor axis, angle phi from A
     *
     * @param A_new
     * @param center
     * @param id
     */
    CEllipse2(STGM::CMatrix2d &A, STGM::CVector2d &center, int id) :
       m_center(center), m_A(A), m_a(0), m_b(0), m_phi(0), m_id(id),  m_type(10)
    {
          int n = 2, info = 0;

          double B[4];
          B[0] = m_A[0][0];
          B[1] = m_A[1][0];
          B[2] = m_A[0][1];
          B[3] = m_A[1][1];

          /** eigenvalue decomposition */
          double eval[2] = {0,0};
          real_eval(B,&n,eval,&info);
          //Rprintf("B: %f %f %f %f \n", B[0],B[1],B[2],B[3]);

          m_majorAxis[0] = B[0];
          m_majorAxis[1] = B[1];
          m_minorAxis[0] = B[2];
          m_minorAxis[1] = B[3];

          if(info != 0) {
        	  error("Eigenvalue decomposition (LAPACK routine) failed in `ellipse2` constructor.");
          } else {

           	/* phi is angle in intersecting plane
           	 * always relative to ´x´ axis */
        	double cos_phi = B[0];
        	double sin_phi = B[1];

            /** relative to z axis */
			//double cos_phi = B[1];
			//double sin_phi = B[3];

			if( (cos_phi<0 && sin_phi<0) ||
				(cos_phi<0 && sin_phi>=0))
			{
			  m_phi = atan(sin_phi/cos_phi)+M_PI;
			} else if(cos_phi>0 && sin_phi<0) {
			  m_phi = atan(sin_phi/cos_phi)+2*M_PI;
			} else {
			  m_phi = acos(cos_phi);
			}

		  	/* axes lengths */
			m_b = 1.0/sqrt(eval[1]);
			m_a = 1.0/sqrt(eval[0]);
          }

    }

    /* no re-computation of matrix A */
    CEllipse2(STGM::CVector2d &center, STGM::CMatrix2d &A, STGM::CVector2d &major,
    		   STGM::CVector2d &minor, double a,  double b, double phi, int id) :
		 m_center(center),
		 m_A(A),
		 m_a(a),
		 m_b(b),
		 m_phi(phi),
		 m_id(id),
		 m_type(ELLIPSE_2D),
		 m_majorAxis(major),
		 m_minorAxis(minor)

	{}


    CEllipse2(STGM::CVector2d &center, STGM::CVector2d &major, STGM::CVector2d &minor,
               double a,  double b, int id) :
         m_center(center),
		 m_a(a),
		 m_b(b),
		 m_phi(0),
		 m_id(id),
		 m_type(ELLIPSE_2D),
         m_majorAxis(major),
		 m_minorAxis(minor)

    {
    	ComputeMatrix();
    }

    void ComputeMatrix()
    {
          STGM::CMatrix2d B;
          B[0][0] = m_majorAxis[0];
          B[1][0] = m_majorAxis[1];
          B[0][1] = m_minorAxis[0];
          B[1][1] = m_minorAxis[1];

          /** compute A */
          m_A[0][0] = 1.0 / SQR(m_a);
	      m_A[1][1] = 1.0 / SQR(m_b);

	      m_A = m_A * B;
		  B.Transpose();
		  m_A = B * m_A;

		  /* angle in the intersecting plane is relative to ´x´ axis: */
          double cos_phi = m_majorAxis[0];
          double sin_phi = m_majorAxis[1];

          if( (cos_phi<0 && sin_phi<0) ||
        	  (cos_phi<0 && sin_phi>=0)) {
              m_phi = atan(sin_phi/cos_phi)+M_PI;
          } else if(cos_phi>0 && sin_phi<0) {
              m_phi = atan(sin_phi/cos_phi)+2.0*M_PI;
          } else {
        	  m_phi = acos(cos_phi);
          }
	}

    /**
     * @brief Set dx/dt=0 and dy/dt=0 -> reorder for t values
     *     Helper functions for calculation of extreme points of the ellipse
     *
     * @return Extreme points of the ellipse,
     */
    STGM::CPoint2d getValue() {  return STGM::CPoint2d(atan(-m_b*tan(m_phi)/m_a),atan(m_b/(tan(m_phi)*m_a))); }

    STGM::CPoint2d PointOnEllipse(const double t) {
      return STGM::CPoint2d(m_center[0] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi),
          m_center[1] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi));
    }

    /**
     * @brief Minimum/maximum y coordinates of the ellipse
     *        which is not unique and needs to be sorted
     *
     * @return
     */
    STGM::CPoint2d getMinMax_Y() {
      double t = getValue()[1];
      double y1 = m_center[1] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi);
      double y2 = m_center[1] + m_a*cos(t+M_PI)*sin(m_phi)+m_b*sin(t+M_PI)*cos(m_phi);

      return ( y1<y2 ? STGM::CPoint2d(y1,y2) : STGM::CPoint2d(y2,y1)  );
    }
    /**
     * @brief Minimum/maximum x coordinates of the ellipse
     *       which is not unique and needs to be sorted
     *
     * @return
     */
    STGM::CPoint2d getMinMax_X() {
      double t = getValue()[0];
      double x1 = m_center[0] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi);
      double x2 = m_center[0] + m_a*cos(t+M_PI)*cos(m_phi)-m_b*sin(t+M_PI)*sin(m_phi);

      return ( x1<x2 ? STGM::CPoint2d(x1,x2) : STGM::CPoint2d(x2,x1)  );
    }

   /**
    * @brief Minimum/Maximum coordinates to run through
    *        for digitizing the ellipse intersections
    *
    * @return [min_x,max_x],[min_y,max_y]
    */
    PointVector2d getMinMaxPoints() {
      PointVector2d p;
      p.push_back(getMinMax_X());
      p.push_back(getMinMax_Y());
      return p;
    };

    void move(CVector2d &v) {  m_center -= v; }


    bool isInside(double x, double y)  {
     double d1 = cos(m_phi)*(x-m_center[0]) + sin(m_phi)*(y-m_center[1]);
     double d2 = sin(m_phi)*(x-m_center[0]) - cos(m_phi)*(y-m_center[1]);

     return (pow(d1 ,2.0)/pow(m_a,2) + pow(d2 ,2.0)/pow(m_b,2)) <= 1;
    }

    /**
     * @brief   sample the ellipse border with np points
     *
     * @param P         Point vector of sampled points
     * @param np        number of points to sample
     */
    void samplePoints(STGM::PointVector2d &P, int np) {
      double t=0.0, s=2.0*M_PI / (double)np;
      for(int k=0; k<np; k++) {
          P.push_back( PointOnEllipse(t) );
          t += s;
      }
    }


    /**
     *
     * @return major axis
     */
    double a() const { return m_a; }
    /**
     *
     * @return minor axis
     */
    double b() const { return m_b; }

    /**
     *
     * @return angle major axis to x axis
     */
    double phi() const { return m_phi; }
    double area() const { return M_PI * m_a *m_b; }

    STGM::CVector2d &center() { return m_center;}
    const STGM::CVector2d &center()  const { return m_center;}

    /**
     * @return Ellipse matrix
     */
    const CMatrix2d &MatrixA() const { return m_A; }

    inline int Id() { return m_id; }

    /**
     * Ellipse2d type
     * @return
     */
    inline int getType() { return m_type; }

    /**
    *
    * @return bounding rectangle
    */
    CBoundingRectangle & getBoundingRectangle() { return m_br; }

    /**
    * @brief Does the ellipse intersect the window or
    *        is it fully contained ?
    *
    * @param win Window
    * @return
    */

    /**
    * @brief Extreme points on the ellipse
    *        Not unique, Maximum can return minimum coordinate point
    * @return Point
    */
    STGM::CVector2d getMaxEllipsePoint_X() {
      return STGM::CVector2d(getMinMax_X()[0],getMinMax_Y()[0]);
    }

    STGM::CVector2d getMaxEllipsePoint_Y() {
      return STGM::CVector2d(getMinMax_X()[0],getMinMax_Y()[1]);
    }

    STGM::CVector2d getMinEllipsePoint_X() {
      return STGM::CVector2d(getMinMax_X()[1],getMinMax_Y()[0]);
    }

    STGM::CVector2d getMinEllipsePoint_Y() {
      return STGM::CVector2d(getMinMax_X()[1],getMinMax_Y()[1]);
    }

    /** only test if center is within window */
    inline bool isInWindow(STGM::CWindow &win) {
      if( win.PointInWindow(m_center) == 0)
       return true;
     return false;
    }

    STGM::CVector2d &majorAxis()  { return m_majorAxis; }
    const STGM::CVector2d &majorAxis()  const { return m_majorAxis;}

    STGM::CVector2d &minorAxis()  { return m_minorAxis;}
    const STGM::CVector2d &minorAxis()  const { return m_minorAxis;}

  private:
    STGM::CVector2d m_center;
    CMatrix2d m_A;
    double m_a, m_b, m_phi;
    int m_id, m_type;
    CBoundingRectangle m_br;
    STGM::CVector2d m_majorAxis, m_minorAxis;
  };


  /**
   * @brief Ellipse in 3d
   */

  class CEllipse3 : public CGeometry
  {
  public:

      CEllipse3() :
        m_center(STGM::CVector3d(0,0,0)),
		m_n(STGM::CVector3d(0,0,1)),
        m_majorAxis(STGM::CVector3d(0,0,1)),
		m_minorAxis(STGM::CVector3d(0,0,0)),
		m_plane(STGM::CPlane()),
        m_a(1),
		m_b(1),
		m_phi(0),
		m_i(0), m_j(1),
		m_type(7),   // default is a full ellipse in 3D
		m_side(0),
		m_side0(0)   // default values for side of caps relative to the plane (not intersected yet)
      {
         m_psi[0] = 0;
         m_psi[1] = 0;
         setPlaneIdx();
      }

      virtual ~CEllipse3() {};

      CEllipse3(STGM::CVector3d &center, STGM::CVector3d &n,
                STGM::CVector3d &major,STGM::CVector3d &minor,
				STGM::CVector3d &mPoint0, STGM::CVector3d &mPoint1,
                double a, double b, double phi, double psi0, double psi1,
				double r0, double r1, int type, int side = 0) :
         m_center(center),
		 m_n(n),
		 m_majorAxis(major),
		 m_minorAxis(minor),
		 m_plane(STGM::CPlane(n)),
		 m_a(a), m_b(b),
         m_phi(phi), m_i(0), m_j(1),
         m_type(type),
		 m_side(side),
		 m_side0(0)									     	// default values because not intersected yet
      {
        m_psi[0] = psi0;
        m_psi[1] = psi1;
        setPlaneIdx();										// set m_i and m_j to select coordinates

        if(m_type == ELLIPSE_SEGMENT)
        {
        	m_circle1 = CCircle3(mPoint0,r0,n,1);
        	m_circle2 = CCircle3(mPoint1,r1,n,1);
        	setReferenceSide();								// set m_side0

        } else if(m_type == ELLIPSE_ARC || m_type == CAP)
        {
        	m_circle1 = CCircle3(mPoint0,r0,n,1);
        	setReferenceSide();
        }

      };


      inline void setPlaneIdx() {
		  switch(m_plane.idx()) {
			case 0: m_i=1; m_j=2; break; // YZ
			case 1: m_i=0; m_j=2; break; // XZ
			case 2: m_i=0; m_j=1; break; // XY
		  }
	  }

      /**
       *
       * @return major axis
       */
      double a() const { return m_a; }
      /**
       *
       * @return minor axis
       */
      double b() const { return m_b; }

      /**
       *
       * @return angle major axis to x axis
       */
      double & phi() { return m_phi; }
      double phi() const { return m_phi; }

      /**
       *
       * @return
       */
      const double *psi() const { return m_psi; }

      STGM::CVector3d &center()  { return m_center;}
      const STGM::CVector3d &center()  const { return m_center;}

      STGM::CVector3d &majorAxis()  { return m_majorAxis;}
      const STGM::CVector3d &majorAxis()  const { return m_majorAxis;}

      STGM::CVector3d &minorAxis()  { return m_minorAxis;}
      const STGM::CVector3d &minorAxis()  const { return m_minorAxis;}


      /**
       * @brief Set dx/dt=0 and dy/dt=0 -> reorder for t values
       *     Helper functions for calculation of extreme points of the ellipse
       *
       * @return Extreme points of the ellipse
       */
      STGM::CPoint2d getValue() {
        return STGM::CPoint2d(atan(-m_b*tan(m_phi)/m_a),atan(m_b/(tan(m_phi)*m_a)));
      }

      STGM::CPoint2d PointOnEllipse(const double t) {
        return STGM::CPoint2d(m_center[m_i] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi),
                              m_center[m_j] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi));
      }

      /**
      * @brief Minimum/maximum y coordinates of the ellipse
      *        which is not unique and needs to be sorted
      *
      * @return
      */
     STGM::CPoint2d getMinMax_Y() {
       double y1,y2,t;
       t = getValue()[1];
       y1 = m_center[m_j] + m_a*cos(t)*sin(m_phi)+m_b*sin(t)*cos(m_phi);
       y2 = m_center[m_j] + m_a*cos(t+M_PI)*sin(m_phi)+m_b*sin(t+M_PI)*cos(m_phi);

       return ( y1<y2 ? STGM::CPoint2d(y1,y2) : STGM::CPoint2d(y2,y1)  );
     }
     /**
      * @brief Minimum/maximum x coordinates of the ellipse
      *       which is not unique and needs to be sorted
      *
      * @return
      */
     STGM::CPoint2d getMinMax_X()
     {
        double x1,x2,t;
        t = getValue()[0];
        x1 = m_center[m_i] + m_a*cos(t)*cos(m_phi)-m_b*sin(t)*sin(m_phi);
        x2 = m_center[m_i] + m_a*cos(t+M_PI)*cos(m_phi)-m_b*sin(t+M_PI)*sin(m_phi);

        return ( x1<x2 ? STGM::CPoint2d(x1,x2) : STGM::CPoint2d(x2,x1)  );
     }

     void move(CVector2d &v) {
    	 m_center[m_i] -= v[0];
    	 m_center[m_j] -= v[1];
    	 if(m_type > 7) {
    		 m_circle1.move(v);
    		 m_circle2.move(v);
    	 }
     }

     bool isInsideEllipse(double x, double y) {
       double d1 = cos(m_phi)*(x-m_center[m_i]) + sin(m_phi)*(y-m_center[m_j]);
       double d2 = sin(m_phi)*(x-m_center[m_i]) - cos(m_phi)*(y-m_center[m_j]);
       return (pow(d1 ,2.0)/pow(m_a,2) + pow(d2 ,2.0)/pow(m_b,2)) <= 1;
     }

    /**
     * @brief Is requested point inside the ellipse / ellipse_arc / ellipse_segment ?
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @return boolean
     */
     bool isInside(double x, double y);

     /**
      * @brief Does the ellipse/ellipse arc/ellipse segment fully lie inside the window?
      *
      * @param win
      * @return true/false
      */
     bool isInWindow(STGM::CWindow &win);

     /**
      * @brief Minimum/Maximum coordinates to run through
      *        for digitizing the ellipse intersections
      *
      * @return [min_x,max_x],[min_y,max_y]
      */
     PointVector2d getMinMaxPoints();

     /**
      * @brief Minimum/Maximum coordinates of points on the ellipse
      *
      * @return All four points, not sorted
      */
     PointVector2d getEllipseExtremePoints();

     /**
      * @brief Extreme points on the ellipse
      *        Not unique, Maximum can return minimum coordinate point
      * @return Point
      */
     STGM::CPoint2d getMaxEllipsePoint_X() {
      double t = getValue()[0];
      return PointOnEllipse(t);
     }

     STGM::CPoint2d getMaxEllipsePoint_Y() {
      double t = getValue()[1];
      return PointOnEllipse(t);
     }

     STGM::CPoint2d getMinEllipsePoint_X() {
      double t = getValue()[0]+M_PI;
      return PointOnEllipse(t);
     }

     STGM::CPoint2d getMinEllipsePoint_Y() {
      double t = getValue()[1]+M_PI;
      return PointOnEllipse(t);
     }

     /**@brief Get the side of the cut off vector on which the point lies
      *
      * @param p Point to determine the side of
      * @param idx Either idx=0 for first cut off angle (ELLIPSE_ARC) and idx=1 for the second (ELLIPSE_SEG)
      * @return -1 or 0  or 1
      */
     int whichSide(STGM::CPoint2d p, int idx ) {
       STGM::CPoint2d p1(PointOnEllipse(m_psi[idx]));
       STGM::CPoint2d p2(PointOnEllipse(2*M_PI-m_psi[idx]));

       return SGN( m_minorAxis[m_j]*(p[0]-p2[0])-m_minorAxis[m_i]*(p[1]-p2[1]));
     }

     /**
      * @brief Set the side, where extreme points
      *        of the ellipse are not cut off by
      *        the angle psi[0]. This function is
      *        only once called by the intersector
      *        for ellipse type {ELLIPSE_ARC | ELLIPSE_SEG }.
      *
      */
     void setReferenceSide() {
       // choose either 0
       double t = 0;
       // or PI as a point on the ellipse to determine the the cut off side
       if(m_side<0) t = M_PI;
       m_side0 = whichSide(PointOnEllipse(t),0);
     }

     /**
      *
      * @return bounding rectangle
      */
     CBoundingRectangle & getBoundingRectangle() { return m_br; }

       /**
      * @brief Comparison functions
      */
     static bool compareX (STGM::CPoint2d x,STGM::CPoint2d y) { return (x[0]<y[0]);  }
     static bool compareY (STGM::CPoint2d x,STGM::CPoint2d y) { return (x[1]<y[1]);  }


     /** members */
     STGM::CVector3d m_center, m_n;
     STGM::CVector3d m_majorAxis, m_minorAxis;
     STGM::CPlane m_plane;

     double m_a, m_b;
     double m_phi, m_psi[2];
     int m_i, m_j;
     int m_type;                                // ellipse type
     int m_side;                                // side of origin0 to the plane
     int m_side0;                               // side, where extreme points of the ellipse are not cut off by angle psi
     CCircle3 m_circle1, m_circle2;       // circle caps if type {ELLIPSE_ARC | ELLIPSE_SEG }
     CBoundingRectangle m_br;
  };

  /**
   * @brief disc projection  for spheroid and cylinder
   *
   * @param center     3d center of spheroid or clyinder
   * @param u          orientation vector
   * @param a          axis length (spheroid) or radius (cylinder)
   * @param phi        angle phi
   * @param id         id
   */
  STGM::CEllipse2 crackProjection(STGM::CVector3d &center, STGM::CVector3d &u, double a, double phi, int id);

   /**
    *   @brief CSpheroid class
    */
   class CSpheroid : public CGeometry{
     public:
      const char *m_label;

      CSpheroid(CVector3d center, double a, double c, double b, CVector3d u,
                 double theta, double phi, int id, const char *label="N", int interior=1)
      : m_label(label),
        m_center(center),
        m_u(u),
        m_a(a), m_b(b), m_c(c), 		/* a,c are 1st and 2nd semi-minor lengths, b is semi-major length (conform with FBA implementation) */
	    m_theta(theta),
        m_phi(phi),
        m_id(id),
        m_crack(0),
        m_interior(interior)
      {
        m_R = RotationMatrixFrom001(u);
        m_u.Normalize();
        ComputeMatrixA();
      }

       /**
        * @return rotation matrix
        */
       const CMatrix3d & rotationMatrix() const { return m_R; }

       double distance(CSpheroid &s) {  return spheroidDistanceAsCylinder(s);  }

       /**
        * @brief Approximate (euclidian) distance of spheroids
        *        by min distance of cylinders
        *
        * @param sp  spheroids
        * @return
        */
       double spheroidDistanceAsCylinder(CSpheroid &sp) const;

       /**
        *
        * @return the crack projection
        */
       CEllipse2 delamProjection() const;

       /**
       * @return orthogonal projection of spheriod
       */
       CEllipse2 spheroidProjection() const;
      /**
       *
       * @return orthogonal projection of spheriod's minor axis circle
       */
       CEllipse2 crackProjection() const;

       /**
        * @brief Sample points and return area of projection
        *
        * @param P       Point vector of points
        * @param np      number of points to sample
        * @return        area of ellipse
        */
       inline double projectedPointsWithArea(PointVector2d &P, int np) {
           CEllipse2 e = spheroidProjection();
           e.samplePoints(P,np);
           return e.area();
       }

       /**
        * @brief The defining spheroid 3d matrix
        */
       void ComputeMatrixA();
       const CMatrix3d &MatrixA() const { return m_A; }

       CVector3d &center() { return m_center; }
       const CVector3d &center() const { return m_center; }

       double a() const { return m_a; }
       double c() const { return m_c;  }
       double b() const { return m_b;  }

       CVector3d& u() { return m_u; }
       const CVector3d& u() const { return m_u; }

       double phi()   const { return m_phi; }
       double theta() const { return m_theta; }

       inline int Id() { return m_id; }
       void setCrackType(int crack) { m_crack=crack; }

       int &interior() { return m_interior; }
       int interior() const { return m_interior; }

       const char * label() const { return m_label; }


    private:
      CVector3d m_center, m_u;
      double m_a, m_b, m_c;
      double m_theta, m_phi;
      int m_id, m_crack, m_interior;
      CMatrix3d m_R, m_A, m_invA;

   };

   /**
    *  @brief Cylinder
    */

    class CCylinder : public CGeometry
    {
    public:
      const char *m_label;

      CCylinder(CVector3d &center, CVector3d &u, double h, double r,
                 double theta, double phi, int id, const char *label="N", int interior=1 ) :
            m_label(label),
            m_center(center),
            m_u(u),
            m_h(h),
            m_r(r),
            m_theta(theta),
            m_phi(phi),
            m_id(id),
            m_interior(interior),
            m_crack(0)
        {
          m_R = RotationMatrixFrom001(u);
          m_u.Normalize();
          updateOrigins();
        };

        ~CCylinder() {};

        inline int Id() { return m_id; }

        /**
        * @return Center of the Cylinder.
        */
        CVector3d & center() { return m_center; }
        const CVector3d & center() const { return m_center; }

        double cylinderDistance(CCylinder &cyl) const;
        inline double distance(CCylinder &s) {  return cylinderDistance(s); }

        /**
        * @return Origin0 of the Cylinder.
        */
        CVector3d &origin0() { return m_origin0; }
        const CVector3d &origin0() const { return m_origin0; }

        /**
        * @return Origin1 of the Cylinder.
        */
        CVector3d &origin1() { return m_origin1; }
        const CVector3d &origin1() const { return m_origin1; }

        /**
        * @return rotation matrix
        */
        const CMatrix3d &rotationMatrix() const { return m_R; }

        /**
        * @return Rotation axis direction vector.
        */
        CVector3d &u() { return m_u; }
        const CVector3d &u() const { return m_u; }

        /**
        * @return radius of cylinder
        */
        double r() const { return m_r; }
        double &r() { return m_r; }

        /**
        * @return Height of the Cylinder.
        */
        double h() const { return m_h; }
        double &h() { return m_h; }

        /**
        * @return Polar angle of the Cylinder.
        */
        double theta() const { return m_theta; }
        double &theta() { return m_theta; }

        /**
        * @return Plane angle of the Cylinder.
        */
        double phi() const { return m_phi; }
        double &phi() { return m_phi; }

        void setCrackType(int crack) { m_crack = crack; }

        CEllipse2 crackProjection() const;

        double delamProjection(PointVector2d &P, int npoints);
        double projectedPointsWithArea(PointVector2d &P, int npoints);

        int interior() const { return m_interior; }
        int &interior() { return m_interior; }

        const char * label() const { return m_label; }

        /**
         * @brief Update origins of cylinder after rotation
         */
        void updateOrigins() {
          double len = 0.5*m_h;
          m_origin0[0] = m_center[0]-len*m_u[0];
          m_origin0[1] = m_center[1]-len*m_u[1];
          m_origin0[2] = m_center[2]-len*m_u[2];

          m_origin1[0] = m_center[0]+len*m_u[0];
          m_origin1[1] = m_center[1]+len*m_u[1];
          m_origin1[2] = m_center[2]+len*m_u[2];
        }


    private:
        // members
        CVector3d m_center, m_u;
        CVector3d m_origin0, m_origin1;
        CMatrix3d m_R;
        double m_h, m_r;
        double m_theta, m_phi;
        int m_id, m_interior, m_crack;
    };


   /**
     *  @brief Box
     */
    class CBox3 : public CGeometry
    {
    public:

      CBox3 ()
       : m_center(STGM::CVector3d(0.5,0.5,0.5)),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(STGM::CPoint3d(1,1,1)),
         m_low(0,0,0),m_up(1,1,1)
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5;
        m_extent[1] = 0.5;
        m_extent[2] = 0.5;

        ConstructBoundingPlanes();
      }

      CBox3 (double xrange[2], double yrange[2], double zrange[2])
       :  m_u(STGM::CVector3d(1.0,0.0,0.0)),
          m_v(STGM::CVector3d(0.0,1.0,0.0)),
          m_w(STGM::CVector3d(0.0,0.0,1.0)),
          m_size(STGM::CPoint3d(fabs(xrange[1]-xrange[0]),fabs(yrange[1]-yrange[0]),fabs(zrange[1]-zrange[0]))),
          m_low(xrange[0],yrange[0],zrange[0]),
          m_up(xrange[1],yrange[1],zrange[1])
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5*m_size[0];
        m_extent[1] = 0.5*m_size[1];
        m_extent[2] = 0.5*m_size[2];

        m_center[0] = xrange[1]-m_extent[0];
        m_center[1] = yrange[1]-m_extent[1];
        m_center[2] = zrange[1]-m_extent[2];

        ConstructBoundingPlanes();
      }

      CBox3 (STGM::CVector3d center, STGM::CPoint3d size)
       : m_center(center),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(size)
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        m_extent[0] = 0.5*m_size[0];
        m_extent[1] = 0.5*m_size[1];
        m_extent[2] = 0.5*m_size[2];

        m_low[0] = m_center[0]-m_extent[0];
        m_low[1] = m_center[1]-m_extent[1];
        m_low[2] = m_center[2]-m_extent[2];

        m_up[0] = m_center[0]+m_extent[0];
        m_up[1] = m_center[1]+m_extent[1];
        m_up[2] = m_center[2]+m_extent[2];

        ConstructBoundingPlanes();
      }

      /** @brief Constructor Box from left lower point (0,0,0)
       *
       * @param a  extent in x in first direction
       * @param b  extent in y in second direction
       * @param c  extent in y in third direction
       */

      CBox3 (double a, double b, double c)
      :  m_center(STGM::CVector3d(0.5*a,0.5*b,0.5*c)),
         m_u(STGM::CVector3d(1.0,0.0,0.0)),
         m_v(STGM::CVector3d(0.0,1.0,0.0)),
         m_w(STGM::CVector3d(0.0,0.0,1.0)),
         m_size(STGM::CPoint3d(a,b,c)),          /* length in each direction */
         m_low(0,0,0),m_up(a,b,c)
      {
        m_axis[0] = &m_u;
        m_axis[1] = &m_v;
        m_axis[2] = &m_w;

        /* half lengths */
        m_extent[0] = 0.5*a;
        m_extent[1] = 0.5*b;
        m_extent[2] = 0.5*c;

        ConstructBoundingPlanes();
      }

      virtual ~CBox3 () {};

      // public members
      void setExtent(double a, double b, double c);

      double volume() const { return  m_size[0]*m_size[1]*m_size[2]; };

      void ConstructBoundingPlanes();

      void ConstructBoxLateralPlanes();

      /**\brief Construct and return lateral planes
       *
       * @return lateral planes
       */
       LateralPlanes & getLateralPlanes() {
        if(m_lateral_planes.empty())
          ConstructBoxLateralPlanes();
        return m_lateral_planes;
      };

      double PointInBox3(STGM::CVector3d &point);

      const std::vector<CPlane> &getPlanes() const { return m_planes; };

      /**
       * @brief Get planes perpendicular to the argument plane
       *
       * @param plane
       * @return planes
       */
      const std::vector<CPlane> getPlanes( CPlane &plane) const {
        std::vector<CPlane> planes;
        int l = plane.idx();
        for(int k=0; k<4; ++k) {
            if(k!=2*l && k!=2*l+1)
              planes.push_back(m_planes[k]);
        }
           return planes;
      }

      STGM::CVector3d m_center;
      STGM::CVector3d m_u, m_v, m_w, *m_axis[3];
      double m_extent[3];

      STGM::CPoint3d m_size;
      STGM::CVector3d m_low, m_up;
      std::vector<CPlane> m_planes, m_lateral_planes;
    };


    // typedefs 3D
    typedef std::vector<CSphere> Spheres;
    typedef std::vector<CSpheroid> Spheroids;
    typedef std::vector<CCylinder> Cylinders;

    // typedefs 2D
    typedef std::vector<CEllipse2> Ellipses2;
    typedef std::vector<CEllipse3> Ellipses3;


} /* namespace STGM */

#endif /* GEOMETRIC_PRIMITIVES_H_ */
