/*
 * IntersectorSpheroidPlane.h
 *
 *  Created on: 31.07.2014
 *      Author: franke
 */

#ifndef INTERSECTOR_H_
#define INTERSECTOR_H_

#include "GeometricPrimitives.h"

namespace STGM
{

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

  template<class T>
  class Intersector
    {
    public:

        virtual ~Intersector () {};

        virtual bool TestIntersection () = 0;
        virtual bool FindIntersection () = 0;

    protected:
        Intersector () {}
    };


  /**
   * @brief IntersectorSpheroidPlane
   *            Only find first intersection with some plane
   */
  template<>
  class Intersector<STGM::CSpheroid>
  {
   public:
    Intersector(CSpheroid &_spheroid, CPlane &_plane, CPoint3d &_dimensions)
     : m_spheroid(_spheroid), m_plane(_plane), m_size(_dimensions), m_type(EMPTY)
    {
    };

    Intersector(CSpheroid &_spheroid, CPoint3d &_dimensions)
        : m_spheroid(_spheroid), m_size(_dimensions), m_type(EMPTY)
    {
    };

    virtual
    ~Intersector() {}

    bool TestIntersection ();
    bool FindIntersection ();

    void setPlane(CPlane &_plane) { m_plane = _plane; }
    CPlane & getPlane() { return m_plane; }

    CEllipse2 & getEllipse() { return m_ellipse; }
    const CEllipse2 & getEllipse() const { return m_ellipse; }

    CSpheroid &getSpheroid() { return m_spheroid; }
    const CSpheroid &getSpheroid() const { return m_spheroid; }

    CGeometry * getObject() { return &m_ellipse; }
    const CGeometry * getObject() const { return &m_ellipse; }

    /**
     * @brief Only check if Spheroid intersects a given plane
     *        and store the intersecting plane, translate center coordinate
     *        periodically to the opposite plane
     *
     * @param  plane   [IN]
     * @return boolean [OUT]
     */
    bool operator() (const CPlane &_plane) {
      m_plane = _plane;
      if(TestIntersection()) {
        int l = m_plane.idx();
        m_spheroid.center()[l] +=  m_size[l] * m_plane.n[l];
        return true;
      }
      return false;
    }

    int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
      int interior = 1;
      for(size_t j=0; j<planes.size() ; ++j) {
          if( operator()(planes[j])) {
            interior=0;
            break;
          }
      }
      return interior;
    }

  private:
    CSpheroid m_spheroid;
    CPlane m_plane;

    CPoint3d m_size;
    IntersectionType m_type;

    CEllipse2 m_ellipse;

  };

  /**
   * \brief Intersector sphere
   */
  template<>
  class Intersector<STGM::CSphere>
  {
  public:
    Intersector(CSphere& sphere, CPlane& plane, CPoint3d size) :
       m_sphere(sphere), m_plane(plane), m_size(size), m_type(EMPTY)
    {}

    Intersector(CSphere &_sphere, CPoint3d &_dimensions) :
       m_sphere(_sphere), m_size(_dimensions), m_type(EMPTY)
    {
    };

    virtual
    ~Intersector() {}

    bool TestIntersection ();
    bool FindIntersection ();

    void setPlane(CPlane &plane) { m_plane = plane; }
    CPlane & getPlane() { return m_plane; }

    CCircle3 & getCircle() { return m_circle; }
    const CCircle3 & getCircle() const { return m_circle; }

    CSphere &getSphere() { return m_sphere; }
    const CSphere &getSphere() const { return m_sphere; }

    CGeometry * getObject() { return &m_circle;; }
    const CGeometry * getObject() const { return &m_circle; }

    bool operator() (const CPlane &plane) {
      m_plane = plane;
      if(TestIntersection()) {
        int l = m_plane.idx();
        m_sphere.center()[l] +=  m_size[l] * m_plane.n[l];
        return true;
      }
      return false;
    }

    int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
      int interior = 1;
      for(size_t j=0; j<planes.size() ; ++j) {
          if( operator()(planes[j])) {
            interior=0;
            break;
          }
      }
      return interior;
    }

  private:
    CSphere m_sphere;
    CPlane m_plane;

    CPoint3d m_size;
    IntersectionType m_type;

    CCircle3 m_circle;
  };

  template<>
  class Intersector<STGM::CCylinder>
    {
    public:

      Intersector( CCylinder &_cylinder,  CPlane &_plane, CPoint3d &_dimensions)
        : m_cylinder(_cylinder), m_plane(_plane), m_size(_dimensions), m_type(EMPTY), m_side(0)
      {
        setPlaneIdx();
      };


      Intersector(  CCylinder &_cylinder, CPoint3d &_dimensions)
        : m_cylinder(_cylinder), m_plane(STGM::CPlane()), m_size(_dimensions), m_type(EMPTY),  m_side(0)
      {
        setPlaneIdx();
      };


      virtual
      ~Intersector() {};

      CCylinder & getCylinder () { return m_cylinder; };
      const CCylinder & getCylinder () const { return m_cylinder; };

      CCircle3 & getCircle1 () { return m_circle1; };
      const CCircle3 & getCircle1 () const { return m_circle1; };

      CCircle3 & getCircle2 () { return m_circle2; };
      const CCircle3 & getCircle2 () const { return m_circle2; };

      CEllipse3 & getEllipse () { return m_ellipse; };
      const CEllipse3 & getEllipse () const { return m_ellipse; };

      int getType() const { return m_type; };
      int getSide() const { return m_side; };

      CPlane & getPlane() { return m_plane; }

      CGeometry * getObject() {
        if(m_type == CAP || m_type == DISC) {
             return & m_circle1;
         } else {  					/* if(type == ELLIPSE ||  type == ELLIPSE_ARC ||  type == ELLIPSE_SEGMENT) */
             return &m_ellipse;
         }
      }

      /**
       * @param plane
       * @return
       */
      bool TestIntersection(const CPlane &plane);

      /**
       * @return
       */
      bool TestIntersection () {
        return TestIntersection(m_plane);
      }
      /**
       * @return
       */
      bool FindIntersection ();

      /**
       * @return
       */
      IntersectionType FindIntersectionType();

      inline void setPlaneIdx() {
          switch(m_plane.idx()) {
           case 0: m_i=1; m_j=2; break;
           case 1: m_i=0; m_j=2; break;
           case 2: m_i=0; m_j=1; break;
          }
      }

      /**
       *
       * @param spherecenter
       * @param sDist
       * @return
       */
      CCircle3 GetCircle(STGM::CVector3d &center, double sDist);

      /**
       *
       * @param a
       * @param cosTheta
       * @return
       */
      int FindMajorAxisIntersection(const double a, const double cosTheta);
      /**
       *
       * @param spherecenter
       * @param ipt
       * @return
       */
      double GetEllipseSegment(STGM::CVector3d center, const STGM::CVector3d &ipt);


      /**
       * @brief Only check if Spheroid intersects a given plane
       *        and store the intersecting plane, translate center coordinate
       *        periodically to the opposite plane
       *
       * @param  plane   [IN]
       * @return boolean [OUT]
       */
      bool operator() (const CPlane &plane) {
        if(TestIntersection(plane)) {
          int l = plane.idx();
          m_cylinder.center()[l] +=  m_size[l] * plane.n[l];
          m_cylinder.origin0()[l] = m_cylinder.center()[l] - 0.5*m_cylinder.h()*m_cylinder.u()[l];
          m_cylinder.origin1()[l] = m_cylinder.center()[l] + 0.5*m_cylinder.h()*m_cylinder.u()[l];
          return true;
        }
        return false;
      }

      int TestBoxIntersection(std::vector<STGM::CPlane> planes) {
		int interior = 1;
		for(size_t j=0; j<planes.size() ; ++j) {
			if( operator()(planes[j])) {
			  interior=0;
			  break;
			}
		}
		return interior;
     }

     private:
      CCylinder m_cylinder;
      CPlane m_plane;

      CPoint3d m_size;
      IntersectionType m_type;

      int m_side, m_i, m_j;
      CCircle3 m_circle1, m_circle2;
      CEllipse3 m_ellipse;

    public:
       CVector3d ipt0, ipt1;

    };

    /** some type definitions */
    template<class T>
    struct Intersectors { typedef typename std::vector< Intersector<T> > Type;  };


    /**
     * @brief Digitizer class
     */
  	 class CDigitizer
     {

       public:

        CDigitizer(int *w, int nrow, int ncol, double delta) :
          m_w(w), m_nrow(nrow), m_ncol(ncol), m_delta(delta),
          x(STGM::CPoint2d(0,0)), y(STGM::CPoint2d(0,0))
        {
          /** safer: initialize */
          for(int i=0; i < nrow*ncol; i++)
        	m_w[i]=0;

          m_nr = m_nrow-1;
          m_nc = m_ncol-1;
          m_d  = 0.5*m_delta-1e-6;
        }

        virtual ~CDigitizer() {};

        template<typename T>
        void start(typename Intersectors<T>::Type &objects)
        {
           CGeometry *obj;
           PointVector2d p;
           CBoundingRectangle br;

           for(size_t k=0; k<objects.size();k++)
           {
               obj = objects[k].getObject();
			   p = obj->getMinMaxPoints();

			   x=p[0]; y=p[1];
			   br.m_ymin=std::max(0,(int)((y[0]+m_d)/m_delta)); // y-coordinate is related to row number
			   br.m_xmin=std::max(0,(int)((x[0]+m_d)/m_delta)); // x-coordinate is related to col number
			   br.m_ymax=std::min(m_nr,(int)((y[1]-m_d)/m_delta));
			   br.m_xmax=std::min(m_nc,(int)((x[1]-m_d)/m_delta));

			   for(int i=br.m_ymin;i<(br.m_ymax+1);i++)
			   {
				   for(int j=br.m_xmin;j<(br.m_xmax+1);j++)
				   {
					   /** change i and j for column/row major order */
					   if(!m_w[i+j*m_nrow])
						 if(obj->isInside((j+0.5)*m_delta,(i+0.5)*m_delta))
							 m_w[i+j*m_nrow]=1;
				   }
			   }
           }

        }

        // a template operator for intersection objects
        template<typename T> void operator()(T &object);

      private:

        int *m_w, m_nr, m_nc, m_nrow, m_ncol;
        double m_delta, m_d;

        STGM::CPoint2d x, y;

     };


  	 template<typename T>
	 void CDigitizer::operator ()(T &object)
	 {
  		PointVector2d p = object.getMinMaxPoints();
	    x=p[0]; y=p[1];

	    STGM::CBoundingRectangle br;
	    br.m_ymin=std::max(0,(int)((y[0]+m_d)/m_delta)); // y-coordinate is related to row number
	    br.m_xmin=std::max(0,(int)((x[0]+m_d)/m_delta)); // x-coordinate is related to col number
	    br.m_ymax=std::min(m_nr,(int)((y[1]-m_d)/m_delta));
	    br.m_xmax=std::min(m_nc,(int)((x[1]-m_d)/m_delta));

	    for(int i=br.m_ymin;i<(br.m_ymax+1);i++)
	    {
		   for(int j=br.m_xmin;j<(br.m_xmax+1);j++)
		   {
			   /** change i and j for column/row major order */
			   if(!m_w[i+j*m_nrow])
				 if(object.isInside((j+0.5)*m_delta,(i+0.5)*m_delta))
					 m_w[i+j*m_nrow]=1;
		   }
	    }

	 }


  	 /*
  	  * Used for ...
  	  */
     template<class T>
     void digitize(typename STGM::Intersectors<T>::Type &objects, int *w, int *nPix, double delta)
     {
       STGM::CDigitizer digitizer(w,nPix[0],nPix[1],delta);
       digitizer.start(objects);
     }



} /* namespace STGM */

#endif /* INTERSECTOR_H_ */
