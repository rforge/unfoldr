/**
 * @file Vector.h
 * @date 17.02.2014
 * @author: franke
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <R.h>

#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline int SGN(const T a)
{return (a > 0) - (a < 0); }

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

namespace STGM {


/** Simple point class */
template<size_t N>
struct CPoint {
  typedef double value_t;
  value_t p[N];
  size_t n_;

  CPoint(const value_t a_);
  CPoint(const CPoint &x);
  CPoint(value_t x0 = 0.0, value_t x1 = 0.0, value_t x2 = 0.0);
  CPoint(const std::vector<value_t> &x);
  CPoint(value_t * v);

  CPoint& operator= (const CPoint &p_);
  value_t &operator[](size_t i) { return p[i]; }
  const value_t &operator[](size_t i) const { return p[i]; }

  bool operator== (const CPoint &p_) const;
  bool operator <(const CPoint &q) const;

  void Normalize();
  value_t Length() const;

  const value_t* ptr() const { return &p[0]; }
  value_t* ptr()  { return &p[0]; }

};

template <size_t N>
CPoint<N>::CPoint(const value_t a) : n_(N) {
  for (size_t i=0; i<n_; i++) p[i]=a;
}

template <size_t N>
CPoint<N>::CPoint(const CPoint &x) : n_(N) {
  if(x.n_ != n_)
    error("Length error");
  for (size_t i=0; i<n_; i++)
    p[i] =x.p[i];
}

template <size_t N>
CPoint<N>::CPoint(value_t x0, value_t x1, value_t x2) : n_(N) {
    p[0] = x0;
    if (N > 1) p[1] = x1;
    if (N > 2) p[2] = x2;
    if (N > 3) error("Point not implemented for N > 3");
}


template <size_t N>
CPoint<N>::CPoint(const std::vector<value_t> &x) : n_(N) {
  if(x.size() != n_)
    error("Length error");
  Memcpy(p,x.data(),n_);
}

template <size_t N>
CPoint<N>::CPoint(value_t * v) : n_(N) {
  for (size_t i=0; i<n_; i++)
    p[i] = v[i];
}


template <size_t N>
CPoint<N>& CPoint<N>::operator=(const CPoint<N>& x) {
  if(this != &x) {
      for (size_t i=0; i<n_; i++)
        p[i] = x.p[i];
  }
  return *this;
}

template <size_t N>
bool CPoint<N>::operator==(const CPoint<N>& x) const {
  if(this != &x) {
      for (size_t i=0; i<n_; i++)
        if (p[i] != x.p[i])
          return false;
  }
  return true;
}

template <size_t N>
inline typename CPoint<N>::value_t CPoint<N>::Length() const {
  value_t tmp=0;
  for (size_t i=0; i<n_; i++)
    tmp += SQR(p[i]);
   return std::sqrt(tmp);
}

template <size_t N>
inline void CPoint<N>::Normalize() {
   const double invLen = 1.0 / Length();
   for (size_t i=0; i<n_; i++)
     p[i] *= invLen;
}

template <size_t N>
bool CPoint<N>::operator <(const CPoint &q) const {
     return p[0] < q[0] || (p[0] == q[0] && p[1] < q[1]);
}

typedef CPoint<2> CPoint2d;
typedef CPoint<3> CPoint3d;


/** Simple vector class */
template<class T, size_t N>
class CVector {
  typedef T value_t;

  value_t data[N];
  size_t n_;

public:

  CVector(const T a);
  CVector(const CVector &v);
  CVector(T x0=0.0, T x1=0.0, T x2=0.0);
  CVector(const std::vector<T> &v);
  CVector(T* p);

  CVector& operator= (const CVector &v);
  inline T & operator[](const int i);
  inline const T& operator[](const int i) const;

  inline CVector& operator+=(const CVector& v);
  inline CVector& operator-=(const CVector& v);
  inline CVector& operator*=(const value_t& a);

  size_t size() const { return n_; }

  const value_t* ptr() const { return &data[0]; }
  value_t* ptr()  { return &data[0]; }

  void Normalize();
  value_t Length() const;
  value_t dot( const CVector &u) const;

};


template <class T, size_t N>
inline T & CVector<T,N>::operator[](const int i) { return data[i]; }

template <class T, size_t N>
inline const T & CVector<T,N>::operator[](const int i) const { return data[i]; }


typedef CVector<double,2> CVector2d;
typedef CVector<double,3> CVector3d;


template <class T,size_t N>
CVector<T,N>::CVector(const value_t a) : n_(N) {
  for (size_t i=0; i<n_; i++) data[i]=a;
}

template <class T,size_t N>
CVector<T,N>::CVector(const CVector &v) : n_(N) {
  if(v.size()!=n_)
    error("Vector Length error in copy construct");
  for (size_t i=0; i<n_; i++)
    data[i] = v.data[i];
}


template <class T,size_t N>
CVector<T,N>::CVector(value_t x0, value_t x1, value_t x2) : n_(N) {
    data[0] = x0;
    if (N > 1) data[1] = x1;
    if (N > 2) data[2] = x2;
    if (N > 3) error("Point not implemented for N > 3");
}


template <class T,size_t N>
CVector<T,N>::CVector(const std::vector<value_t> &v) : n_(N) {
  if(v.size() != n_)
    error("Vector Length error in vector copy");
  Memcpy(data,v.data(),n_);
}

template <class T,size_t N>
CVector<T,N>::CVector(T * v) : n_(N) {
  for (size_t i=0; i<n_; i++)
    data[i] = v[i];
}


///unary vector class operators
template <class T,size_t N>
CVector<T,N>& CVector<T,N>::operator=(const CVector<T,N>& v) {
  if(this != &v) {
      for (size_t i=0; i<n_; i++)
        data[i] = v.data[i];
  }
  return *this;
}

/// binary vector class operators
template <class T,size_t N>
inline CVector<T,N> operator- (const CVector<T,N> &u) {
  CVector<T,N> v;
  for (size_t i=0; i<v.size(); i++)
      v[i] = -u[i];
  return v;
}

template <class T,size_t N>
inline CVector<T,N>& CVector<T,N>::operator*=(const value_t &a) {
  for (size_t i=0; i<n_; i++)
      data[i] *= a;
  return *this;
}


template <class T,size_t N>
inline CVector<T,N>& CVector<T,N>::operator+=(const CVector<T,N> &u) {
  for (size_t i=0; i<n_; i++)
      data[i] += u.data[i];
  return *this;
}



template <class T,size_t N>
inline CVector<T,N>& CVector<T,N>::operator-=(const CVector<T,N> &u) {
  for (size_t i=0; i<n_; i++)
      data[i] -= u.data[i];
  return *this;
}


template <class T,size_t N>
inline CVector<T,N> operator*(const double a, const CVector<T,N> &v) {
  CVector<T,N> r(v);
  for (size_t i=0; i<v.size(); i++)
      r[i] *= a;
  return r;
}


template <class T,size_t N>
inline CVector<T,N> operator*(const CVector<T,N> &v, const double a ) {
  CVector<T,N> r(v);
  for (size_t i=0; i<v.size(); i++)
      r[i] *= a;
  return r;
}


template <class T,size_t N>
inline typename CVector<T,N>::value_t CVector<T,N>::dot( const CVector &u) const {
  if(u.size() != n_)
    error("Length error");
  value_t tmp=0;
  for (size_t i=0; i<n_; i++)
    tmp += data[i]*u[i];
   return tmp;
}



template <class T,size_t N>
inline typename CVector<T,N>::value_t CVector<T,N>::Length() const {
  value_t tmp=0;
  for (size_t i=0; i<n_; i++)
    tmp += SQR(data[i]);
   return std::sqrt(tmp);
}

template <class T,size_t N>
inline void CVector<T,N>::Normalize() {
   const double invLen = 1.0 / Length();
   for (size_t i=0; i<n_; i++)
     data[i] *= invLen;
}

inline CVector3d cross(const CVector3d &a, const CVector3d &b) {
  if(a.size() != b.size())
    error("Length error");
  return CVector3d(
      a[1]*b[2]-a[2]*b[1],
      a[2]*b[0]-a[0]*b[2],
      a[0]*b[1]-a[1]*b[0]);
}


/** Simple matrix class */
class CMatrix2d {
  double data[2][2];

  public:
  CMatrix2d() {
    std::memset(data, 0, sizeof(data));
  }

  ~CMatrix2d() {}

  CMatrix2d(const double *a) {
    /**  from R column-major to C row-major */
	data[0][0] = a[0];
	data[0][1] = a[1];
	data[1][0] = a[2];
	data[1][1] = a[3];
  }

  double *operator[](int i) { return data[i]; }
  const double *operator[](int i) const { return data[i]; }

  void Transpose() {
      for (int i = 0; i < 1; ++i) {
        for (int j = i + 1; j < 2; ++j) {
          std::swap(data[i][j], data[j][i]);
        }
      }
  }

  void nullify() { std::memset(data, 0, sizeof(data)); }

};

inline const CMatrix2d operator*(const CMatrix2d &a, const CMatrix2d &b) {
    CMatrix2d m;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        m[i][j] = 0.0;
        for (int k = 0; k < 2; ++k) {
          m[i][j] += a[i][k] * b[k][j];
        }
      }
    }
    return m;
}


class CMatrix3d {
  double data[3][3];

  public:

   CMatrix3d() {
     std::memset(data, 0, sizeof(data));
   }

  ~CMatrix3d() {}

  double *operator[](int i) { return data[i]; }
  const double *operator[](int i) const { return data[i]; }

  void Transpose() {
    for (int i = 0; i < 2; ++i) {
      for (int j = i + 1; j < 3; ++j) {
        std::swap(data[i][j], data[j][i]);
      }
    }
  }

  void nullify() { std::memset(data, 0, sizeof(data));  }

};

/**
 * @param a Left matrix.
 * @param b Right matrix.
 * @return Product of these two matrices.
 */
inline const CMatrix3d operator*(const CMatrix3d &a, const CMatrix3d &b) {
  CMatrix3d c;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      c[i][j] = 0.0;
      for (int k = 0; k < 3; ++k) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return c;
}


}

#endif /* VECTOR_H_ */
