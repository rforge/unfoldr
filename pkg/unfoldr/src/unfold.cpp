/**
 * unfold.cpp
 *
 *  Created on: 18.03.2015
 *      Author: M. Baaske
 */

#include <R_ext/Rdynload.h>

#include "SimSphere.h"
#include "SimEllipsoid.h"
#include "SimCylinder.h"

#include "unfold.h"

#ifdef _OPENMP
 #include <omp.h>
 #include <R_ext/MathThreads.h>
 static int nthreads=1;
#endif

#define DIM_SIZE 6
#define DIM_SIZE_HIST 3

#define ZERO_TOL 10e-9
#define PHI_TOL  10e-13

#define DivPI 1.27323954473516276486

/**  elliptic integral second kind */
extern double elleptint(double, double);

using namespace std;

template<int DIM>
struct  CArray {
  typedef double CTYPE ;

  SEXP A;
  CTYPE* array;
  size_t L,M,N,I,J,K;

  CArray(SEXP _A) : A(_A)
  {
    switch(DIM) {
      case 1:
        error(_("CArray(): Array with dimension less than 2 not implemented!"));
        break;
      case 2:
        init(); I=dims[0]; J=dims[1]; break;
      case 3:
        init(); I=dims[0]; J=dims[1]; K=dims[2]; break;
      case 6:
        init();
        I=dims[0]; J=dims[1];  K=dims[2];
        L=dims[3]; M=dims[4];  N=dims[5];
        J1=J*I; K1=K*J1; L1=L*K1; M1=M*L1;
        break;
      default:
        error(_("CArray(): Array with supplied dimension not implemented!"));
    }
    array=ptr();
    sizeN=size();
  }

  virtual ~CArray() {}

  void init() { int *d=getDims(A); for(int i=0;i<DIM;++i)  dims[i]=d[i]; }

  CTYPE & operator()(int i, int j) { return array[j*I+i]; }
  CTYPE & operator()(int i, int j, int k) { return array[k*I*J+j*I+i]; }
  CTYPE & operator()(int i, int j, int k,int l, int m, int n) { return array[n*M1+m*L1+l*K1+k*J1+j*I+i]; }

  CTYPE * ptr() {  return REAL(A); }
  void operator++(int value) { ++value; }

  size_t indx(int i, int j, int k,int l, int m, int n) {
	  return n*M1+m*L1+l*K1+k*J1+j*I+i;
  }

  void clear() { for(size_t i=0;i<sizeN;i++)  array[i]=0.0; }

  size_t size() const {
    size_t sz=1;
    for(int i=0;i<DIM;++i)
      sz*=dims[i];
    return sz;
  }

  size_t dims[DIM]; // I,J,K,L,M,N
  size_t J1,K1,L1,M1, sizeN;
};

typedef CArray<6> CArray6;
typedef CArray<3> CArray3;



/**
 * @brief Binary search index of element x in vector vec
 *
 * @param v vector of bins
 * @param n length of v
 * @param x element
 * @return index: -1 -> element is outside
 */

int findIndex(double *v, int n, double x) {
    if (x<=v[0]) return -1;
    if (x>v[n-1]) return -1;
    int j,i=0,k=n-1;

    while(j=(k-i)/2+i,!(v[j]<x && x<=v[j+1])) {
      if (x<=v[j]) k=j;
      else i = j+1;
    }
    return j;
}

void
extern_findIndex(double *x, double *v, int *n, int *idx) {
    *idx=findIndex(v,*n,*x);
}

struct K_Oblate_s {
  K_Oblate_s(double *size, double *angle, double *shape) :
    m_A(size), m_alpha(angle), m_S(shape),
    m_a(size), m_Theta(angle), m_s(shape),
    t1(0),t2(0),B(0),tmpS(0)
  {};

  ~K_Oblate_s() {};

  double K(double alpha, double S, double Theta, double s)
  {
     double K1=0,K2=0,phi=0,tmp=0,div=0;

     if(alpha>0 && (M_PI_2-Theta <= alpha && alpha <= M_PI_2))   /* always Theta>0 */
     {
        div=fabs(SQR(sin(alpha))-SQR(cos(Theta)));               /* fabs() for rounding errors only */
        tmp = sqrt(div)/(sin(Theta)*sin(alpha));                 /* always Theta>0 */
        if(!R_FINITE(tmp))
           error(_("KOblate_s(): inf/NaN produced"));
        phi = (fabs(tmp-1.0) < ZERO_TOL ? M_PI_2 : asin(tmp) );  /* asin() domain -> rounding errors */

        if(!R_FINITE(phi)) {
          phi=M_PI_2;
          warning(_("K_Oblate_s(): inf produced"));
        }

        K1=elleptint(phi,B);
     }

     if(s<S && S<=tmpS) {
        phi = asin(sqrt(1.0-SQR(s)/SQR(S))/B);
        if(!R_FINITE(phi)) {
            error(_("K_Oblate_s(): NaN produced"));
        }
        K2=elleptint(phi,B);
      } //else if(S<=s) K2=0;
     else if(S>tmpS) K2=elleptint(M_PI_2,B);

     return MIN(K1,K2);
  }

  inline double operator()(int i, int j, int k,int l, int m, int n)  {
    if(m_a[i+1]<m_A[l+1] || m_a[i+1]<m_A[l])
      return 0;

    B=sin(m_Theta[j+1])*sqrt(1.0-SQR(m_s[k+1]));
    tmpS=m_s[k+1]/sqrt(SQR(m_s[k+1])*SQR(sin(m_Theta[j+1]))+SQR(cos(m_Theta[j+1])));
    t1=DivPI*(m_a[i+1]-sqrt(SQR(m_a[i+1])-SQR(m_A[l+1])));
    t2=DivPI*(m_a[i+1]-sqrt(SQR(m_a[i+1])-SQR(m_A[l])));

    double p=
      t1*K(m_alpha[m+1],m_S[n+1],m_Theta[j+1],m_s[k+1]) - t2*K(m_alpha[m],   m_S[n],   m_Theta[j+1], m_s[k+1])
    + t1*K(m_alpha[m],  m_S[n],  m_Theta[j+1],m_s[k+1]) - t2*K(m_alpha[m+1], m_S[n+1], m_Theta[j+1], m_s[k+1])
    - t1*K(m_alpha[m+1],m_S[n],  m_Theta[j+1],m_s[k+1]) + t2*K(m_alpha[m+1], m_S[n],   m_Theta[j+1], m_s[k+1])
    - t1*K(m_alpha[m],  m_S[n+1],m_Theta[j+1],m_s[k+1]) + t2*K(m_alpha[m],   m_S[n+1], m_Theta[j+1], m_s[k+1]);

    double ret= fabs(p)<ZERO_TOL ? 0.0 : p;
    if(!R_FINITE(ret) || ISNAN(ret)) {
       error(_("K_Oolate_s(): operator return value error."));
    }
    return ret;
  }

  /* members */
   const double *m_A,*m_alpha,*m_S;       // histogram bins, planar
   const double *m_a,*m_Theta,*m_s;       // histogram bins, spatial

   double t1,t2,B,tmpS;
};

struct K_Prolate_s {
   K_Prolate_s(double *size, double *angle, double *shape) :
    m_A(size), m_alpha(angle), m_S(shape),
    m_a(size), m_Theta(angle), m_s(shape),
    t1(0),t2(0),Z(0),Zroot(0),M(0),tmpS(0)
   {
   };

   ~K_Prolate_s() {};

   double K(double alpha, double S, double Theta, double s)
   {
       double K1=0,K2=0,phi=0,tmp=0,eint=0,div=0;
       double eint2 = (M>0 ? elleptint(M_PI_2,M) : M_PI_2);  /* M==0 , if Z==1*/

       if(!R_FINITE(eint2) || ISNAN(eint2))
          error("KProlate_s(): elliptic integral error.");

       if(alpha>0) {
            if(alpha>=Theta) {
                eint=eint2;
            } else {
               phi = asin(cot(Theta)*tan(alpha));
               if(!ISNAN(phi)) {
                 eint=elleptint(phi,M);
               } else {
                 warning(_("KPrblate_s(): NaNs produced"));
                 eint=elleptint(M_PI_2,M);
               }
            }
            K1 = Zroot*eint;
       }  /* else K1=0 */

       if(s<S && S<=tmpS) {
           tmp=sqrt(1.0-SQR(s)/SQR(S));
           div=fabs(tmp-M);
           if(!R_FINITE(div))
             error(_("KProlate(): 'inf/Na/NaN' value produced"));

           if(div<ZERO_TOL) { /* tmp==M */
             eint=eint2;      /* phi=M_PI_2=asin(1) */
           } else {
               phi=asin(tmp/M);
               if(ISNAN(phi))
                 error(_("K_Prolates_s(): NaNs produced."));
               eint=elleptint(phi,M);
           }
           div = fabs(Z-SQR(S)/SQR(s));
           if(!R_FINITE(div))
             error(_("KProlate(): 'inf/Na/NaN' value produced"));
           K2=Zroot*eint;
           if(div>ZERO_TOL) {
               div=sqrt(div);
               if(ISNAN(div))
                 error(_("KProlate(): NaN value produced"));
               K2 -= tmp*div;    /* K2=Zroot*eint-tmp*sqrt(Z-MIN(SQR(S)/SQR(s),Z)); */
           }
           if(!R_FINITE(K2) || ISNAN(K2))
              error("K(): K2 value error.");
       }

       if(S>tmpS)
         return MAX(K1,0.0);

       return MAX(K1+K2-Zroot*eint2,0.0);
   }

   inline double operator()(int i, int j, int k,int l, int m, int n)  {
     if(m_a[i+1]<m_A[l+1] || m_a[i+1]<m_A[l])
       return 0;

     t1 = DivPI*(m_a[i+1]-sqrt(SQR(m_a[i+1])-SQR(m_A[l+1])));
     t2 = DivPI*(m_a[i+1]-sqrt(SQR(m_a[i+1])-SQR(m_A[l])));

     tmpS = sqrt(SQR(m_s[k+1])*SQR(cos(m_Theta[j+1]))+SQR(sin(m_Theta[j+1])));
     Z = 1.0+(1/SQR(m_s[k+1])-1.0)*SQR(sin(m_Theta[j+1]));
     M = sqrt((Z-1.0)/Z);
     Zroot = sqrt(Z);

     double p=
      t1*K(m_alpha[m+1],m_S[n+1], m_Theta[j+1],m_s[k+1]) - t2*K(m_alpha[m],   m_S[n],  m_Theta[j+1], m_s[k+1])
    + t1*K(m_alpha[m],  m_S[n],   m_Theta[j+1],m_s[k+1]) - t2*K(m_alpha[m+1], m_S[n+1],m_Theta[j+1], m_s[k+1])
    - t1*K(m_alpha[m+1],m_S[n],   m_Theta[j+1],m_s[k+1]) + t2*K(m_alpha[m+1], m_S[n],  m_Theta[j+1], m_s[k+1])
    - t1*K(m_alpha[m],  m_S[n+1], m_Theta[j+1],m_s[k+1]) + t2*K(m_alpha[m],   m_S[n+1],m_Theta[j+1], m_s[k+1]);

     double ret= fabs(p)<ZERO_TOL ? 0.0 : p;
     if(!R_FINITE(ret) || ISNAN(ret)) {
         error(_("K_Prolate_s(): operator return value error."));
     }

     return ret;
    }

   /* members */
   const double *m_A,*m_alpha,*m_S;       // histogram bins, planar
   const double *m_a,*m_Theta,*m_s;       // histogram bins, spatial

   double t1,t2,Z,Zroot,M,tmpS;
};

/*
template<class STYPE>
void KFunctor(CArray6 &P, STYPE p) {
  size_t i,j,k,l,m,n;
  for(n=0;n<P.N;++n)
     for(m=0;m<P.M;++m)
        for(l=0;l<P.L;++l)
          for(k=0;k<P.K;++k)
            for(j=0;j<P.J;++j)
               for(i=0;i<P.I;++i)
                 P(i,j,k,l,m,n) =  p(i,j,k,l,m,n);
}
*/

template<class STYPE>
void KFunctor(CArray6 &P, STYPE p) {
  size_t l,m,n,pos;
  size_t K1 = P.K*P.J*P.I;
  double *pP = P.ptr();

#ifdef _OPENMP
if(nthreads==1) {
    for(n=0;n<P.N;++n)
       for(m=0;m<P.M;++m)
          for(l=0;l<P.L;++l) {
              pP = &P(0,0,0,l,m,n);
              for(pos=0;pos<K1;++pos)
               pP[pos] = p(pos%P.I,(pos/P.I)%P.J,pos/(P.I*P.J),l,m,n);
          }

} else {
    for(n=0;n<P.N;++n)
         for(m=0;m<P.M;++m)
              for(l=0;l<P.L;++l) {
                   pP = &P(0,0,0,l,m,n);
                   // transform three nested for loops into one
                   // and calculate indeces i,j,k
                   #pragma omp parallel for num_threads(nthreads) \
                    firstprivate(p,l,m,n) private(pos)
                       for(pos=0;pos<K1;++pos)
                         // i=pos%P.I; j=(pos/P.I)%P.J; k=pos/(P.I*P.J);
                         pP[pos] = p(pos%P.I,(pos/P.I)%P.J,pos/(P.I*P.J),l,m,n);
               }

}
#else
  for(n=0;n<P.N;++n)
    for(m=0;m<P.M;++m)
      for(l=0;l<P.L;++l) {
         pP = &P(0,0,0,l,m,n);
         for(pos=0;pos<K1;++pos)
           pP[pos] = p(pos%P.I,(pos/P.I)%P.J,pos/(P.I*P.J),l,m,n);
      }
#endif
}


template<int KTYPE> struct kernelsp {};
// planar sampling in planar section
template<>
struct kernelsp<0> {
  static double kernelfun(double u, double s) {   if(s<0.0) return u;  if(s<=u)  return sqrt(SQR(u)-SQR(s)); return 0.0;  }
};
// linear sampling in planar section
template<>
struct kernelsp<1> {
  static double kernelfun(double u, double s) {  if(s<0.0) return M_PI/4.0*SQR(u); if(s<=u)  return M_PI/4.0*(SQR(u)-SQR(s)); return 0.0;  }
};

template<int KTYPE>
struct KernFunSP {
  KernFunSP(double *bx) : b(bx) {};
  double operator()(int k, int i) {
      return (kernelsp<KTYPE>::kernelfun(b[i+1], b[k])
             -kernelsp<KTYPE>::kernelfun(b[i+1], b[k+1]));
  }
  double *b;
};

SEXP EMS(SEXP R_P, SEXP R_F, SEXP R_cond) {
  size_t i,j,k,l,m,n,la=0;
  size_t nla = (size_t) asInteger(getListElement( R_cond, "maxSteps"));

#ifdef _OPENMP
 nthreads=MAX(asInteger(getListElement( R_cond, "nCores")),1);
#endif

  if(!isArray(R_P) || length(getAttrib(R_P,R_DimSymbol))!=DIM_SIZE)
      error(_("EMS(): Expected coefficient matrix of type 'array'!"));

  if(!isArray(R_F) || length(getAttrib(R_F,R_DimSymbol))!=DIM_SIZE_HIST)
    error(_("EMS(): Expected the input histogram of type 'array'!"));

  CArray6 P(R_P);                           /* coefficient array */
  SEXP R_H;                                 /* final estimated histogram,  set to input histogram first */

  /** @todo: different bin sizes for spatial histogram */
  PROTECT(R_H=duplicate(R_F));
  CArray3 H(R_H), F(R_F);

  SEXP R_ts, R_rs;
  PROTECT(R_ts=alloc3DArray(REALSXP,F.I,F.J,F.K));
  PROTECT(R_rs=alloc3DArray(REALSXP,F.I,F.J,F.K));

  CArray3 ts(R_ts), rs(R_rs);               /* 3d arrays for t_ijk_sum and r_lmn_sum  */
  double s=0,rsum=0,tsum=0;                 /* temporary sums     */

  /* calculate t_ijk_sum */
  double *p = P.ptr();
  for(k=0; k<F.K; ++k)
    for(j=0; j<F.J; ++j)
      for(i=0; i<F.I; ++i) {
        tsum=0.0;
        p = &P(i,j,k,0,0,0);
        for(n=0; n<F.K; ++n)
          for(m=0; m<F.J; ++m)
            for(l=0; l<F.I; ++l)
              tsum += *(p+n*P.M1+m*P.L1+l*P.K1);         /*  tsum += P(i,j,k,l,m,n);  */
        ts(i,j,k)=tsum;
      }


  /* Main EM iteration */
  size_t K1 = F.K*F.J*F.I;
  double *h = H.ptr(), *f=F.ptr(), *r=rs.ptr();          /* pointer to arrays   */

  for(la=0;la<nla;la++)
  {
      /* calculate r_lmn_sum */
      for(n=0;n<F.K; ++n)
        for(m=0; m<F.J; ++m)
          for(l=0; l<F.I; ++l) {
              rsum = 0.0;
              p = &P(0,0,0,l,m,n);
              h = H.ptr();
#ifdef _OPENMP
              if(nthreads==1) {
                  // sum over i,j,k
                  // originally three nested for loops
                  for(k=0;k<K1;++k)
                    rsum += *p++ * *h++;                /* rs(l,m,n) = P(i,j,k,l,m,n)*H(i,j,k); */
              } else {
               #pragma omp parallel for num_threads(nthreads) \
                  private(k) reduction(+:rsum)
                  for(k=0;k<K1;++k)
                    rsum += p[k] * h[k];
              }
#else
              for(k=0;k<K1;++k)
                rsum += *p++ * *h++;
#endif
             rs(l,m,n) = rsum;
          }

      /** update */
      for(k=0;k<H.K;++k)
        for(j=0;j<H.J;++j)
          for(i=0;i<H.I;++i)
          {
              if(H(i,j,k)<ZERO_TOL)
                 continue;
              r = rs.ptr();
              f = F.ptr();
              p = &P(i,j,k,0,0,0);

              for(n=0,s=0.0;n<H.K; ++n)
                for(m=0; m<H.J; ++m)
                   for(l=0; l<H.I; ++l,r++,f++)
                     /* if(rs(l,m,n)>0.0) s+=P(i,j,k,l,m,n)*F(l,m,n)/rs(l,m,n); */
                     if( *r>0.0) s+=*(p+n*P.M1+m*P.L1+l*P.K1) * *f/ *r;

              if(ts(i,j,k)>0.0)
                H(i,j,k) *= s/ts(i,j,k);
          }
  } /** end loop EM */

  UNPROTECT(3);
  return R_H;
}

void em_saltykov_p(int *nn, double *bin, double *xp) {
  int i,k, n = *nn;
  KernFunSP<0> p(bin); /* up to now: only planar sampling is available */

  for(i=0; i<n;i++)
    for(k=0; k<n;k++)
      *xp++ = p(k,i);  /* P(k,i)=p(k,i); */
}

void em_saltykov(int *_n, int *_nla, double *p, double *y, double *theta)
{
  int i,k,la, n=*_n, nla=*_nla;
  double s, *q=Calloc(n,double), *r=Calloc(n,double);

  for(i=0;i<n;i++)
    for(k=0,q[i]=0.0;k<n;k++) q[i]+=p[i*n+k];

  for(la=0;la<nla;la++)
  {
      for(k=0;k<n;k++)
        for(i=0,r[k]=0.0;i<n;i++)
          r[k]+=p[i*n+k]*theta[i];

      for(i=0;i<n;i++) {
          for(k=0,s=0.0;k<n;k++)
            if(r[k]>0.0) s+=p[i*n+k]*y[k]/r[k];
          if(q[i]>0.0) theta[i]*=s/q[i];
      }
  }

  Free(q);
  Free(r);
}


/**
 * \brief               Calculate coefficient array for EM iteration
 *
 * @param R_A           size, planar
 * @param R_alpha       orientation, planar
 * @param R_S           shape,planar
 * @param R_a           size, spatial
 * @param R_Theta       orientation, spatial
 * @param R_s           shape,spatial
 * @param R_args        type of spheroid ['prolate' | 'oblate']
 *
 * @return              Coefficient array
 */
SEXP CoefficientMatrixSpheroids(SEXP R_A, SEXP R_alpha, SEXP R_S,
                                SEXP R_a, SEXP R_Theta, SEXP R_s, SEXP R_args ) {

#ifdef _OPENMP
 nthreads=MAX(asInteger(VECTOR_ELT(R_args,1)),1);
#endif


  /* spatial class limits equal planar ones in this current implementation */
  SEXP indx;
  PROTECT(indx = NEW_INTEGER(6));
  /* spatial parameters for size, orientation, shape */
  INTEGER(indx)[0] = LENGTH(R_a)-1;  INTEGER(indx)[1] = LENGTH(R_Theta)-1;  INTEGER(indx)[2] = LENGTH(R_s)-1;
  /* and planar parameters */
  INTEGER(indx)[3] = LENGTH(R_A)-1;  INTEGER(indx)[4] = LENGTH(R_alpha)-1;  INTEGER(indx)[5] = LENGTH(R_S)-1;

  SEXP R_P;
  PROTECT(R_P=allocArray(REALSXP,indx));

  CArray6 P(R_P);
  P.clear();

  double *A=REAL(R_A), *S=REAL(R_S), *alpha=REAL(R_alpha);
  /* double *a=REAL(R_a), *s=REAL(R_s), *Theta=REAL(R_Theta); */

  if ( !std::strcmp( translateChar(STRING_ELT(VECTOR_ELT(R_args,0),0)), "prolate" )) {
    K_Prolate_s p(A,alpha,S);
    KFunctor<K_Prolate_s>(P,p);
  } else {
    K_Oblate_s p(A,alpha,S);
    KFunctor<K_Oblate_s>(P,p);
  }

  // set class attributes
  SEXP RClassName = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(RClassName, 0, mkChar("array"));
  SET_STRING_ELT(RClassName, 1, mkChar("PMatrix"));
  classgets(R_P, RClassName);

  UNPROTECT(3);
  return R_P;
}

/**
 *  \brief Binning of vectors,
 *         If some value is lower than lowest bin value findInterval returns zero.
 *         If value is inside some bin: i,j,k=1,...,r with r bin classes
 *
 * @param Rx            vector of values
 * @param Ry            vector of values
 * @param Rz            vector of values
 * @param Rbin_x        bin vector of x
 * @param Rbin_y        bin vector of y
 * @param Rbin_z        bin vector of z
 * @return              Matrix of counts
 */
SEXP Binning3d(SEXP Rx, SEXP Ry, SEXP Rz, SEXP Rbin_x, SEXP Rbin_y, SEXP Rbin_z) {
  int i=0, j=0, k=0;

  /* length of bins */
  int nx=LENGTH(Rx), nxt=LENGTH(Rbin_x), nyt=LENGTH(Rbin_y), nzt=LENGTH(Rbin_z);
  double *x=REAL(Rbin_x), *y=REAL(Rbin_y), *z=REAL(Rbin_z);

  SEXP R_A;
  PROTECT(R_A=alloc3DArray(REALSXP,nxt-1,nyt-1,nzt-1));
  CArray3 A(R_A);
  A.clear();
  for(int p=0; p<nx; ++p) {
    i = findIndex(x,nxt,REAL(Rx)[p]);
    j = findIndex(y,nyt,REAL(Ry)[p]);
    k = findIndex(z,nzt,REAL(Rz)[p]);

    if(i>=0 && j>=0 && k>=0)
      ++A(i,j,k);
  }

#if 0
  for(int k=0;k<A.K;++k) {
       std::cout << ",,"<< k << std::endl;
       for(int i=0;i<A.I;++i) {
         for(int j=0;j<A.J;++j) {
             std::cout << A(i,j,k) << " ";
         }
         std::cout << std::endl;
        }
        std::cout << std::endl;
  }
#endif

  SEXP RClassName = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(RClassName, 0, mkChar("array"));
  SET_STRING_ELT(RClassName, 1, mkChar("triHist"));
  classgets(R_A, RClassName);

  UNPROTECT(2);
  return R_A;
}

/**
 * \brief 1d binning
 *
 * @param Rx   vector of numeric values
 * @param Rbin bin vector
 * @return R numeric vector of count data
 */
SEXP Binning1d(SEXP Rx, SEXP Rbin) {
  int i,k, nt=LENGTH(Rbin);
  SEXP R_A;
  PROTECT(R_A = allocVector(REALSXP,nt-1));

  double *v=REAL(Rbin), *a=REAL(R_A), *x=REAL(Rx);

  for(k=0;k<nt-1;k++) a[k]=0;
  for(k=0;k<LENGTH(Rx);k++) {
      i = findIndex(v,nt,x[k]);
        if(i>=0 )  ++a[i];
  }
  UNPROTECT(1);
  return R_A;
}


/* R Interface functions  */
#define CALLDEF(name, n)  { #name, (DL_FUNC) &name, n}

R_NativePrimitiveArgType myC_t[] = { INTSXP, INTSXP, REALSXP, REALSXP, REALSXP };
R_NativePrimitiveArgType myC_t2[] = { INTSXP, REALSXP, REALSXP };


static R_CMethodDef CEntries[]  = {
	{"em_saltykov", (DL_FUNC) &em_saltykov, 5, myC_t},
	{"em_saltykov_p", (DL_FUNC) &em_saltykov_p, 3, myC_t2},
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef CallEntries[] = {
      CALLDEF(EllipsoidSystem,2),
      CALLDEF(CylinderSystem,2),
      CALLDEF(IntersectSpheroidSystem,6),
	  CALLDEF(IntersectCylinderSystem,6),
      CALLDEF(IntersectSphereSystem,6),
      CALLDEF(UpdateIntersections,2),
      CALLDEF(SphereSystem,2),
      CALLDEF(SimulateSpheresAndIntersect,2),
      CALLDEF(SimulateSpheroidsAndIntersect,2),
	  CALLDEF(SimulateCylindersAndIntersect,2),
      CALLDEF(DigitizeProfiles,3),
	  CALLDEF(Binning3d,6),
      CALLDEF(Binning1d,2),
      CALLDEF(CoefficientMatrixSpheroids,7),
      CALLDEF(EMS,3),
      {NULL, NULL, 0}
};


void R_init_unfoldr(DllInfo *info) {
  R_registerRoutines(info, CEntries,CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

/*
void R_unload_unfold(DllInfo *info){
  // Release resources
}
*/
