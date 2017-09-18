/**
 * @file ellipt.cpp
 * @date 2015-06-1
 *
 * @brief Functions to calculate the complete elliptic
 *        integral of the second kind:
 *        The implementation is based on Carlson's algorithms
 *        and the duplication theorem.
 *
 * @author: M. Baaske
 */

#include "R.h"
#include "Rmath.h"
#include "Vector.h"

using namespace std;

static double D1=3.0/14.0, D2=1.0/6.0, D3=9.0/22.0,
              D4=3.0/26.0, D5=0.25*D3, D6=1.5*D4;

static double ERRTOL_RD=0.001;
static double SMALL_RD=2.0*R_pow(DBL_MAX,-2./3.);
static double BIG_RD=0.1*ERRTOL_RD*R_pow(DBL_MIN,-2./3.);

static double ERRTOL_RF=0.002,
              ONE_THIRD=1.0/3.0,
              C1=1.0/24.0,C2=0.1,C3=3.0/44.0, C4=1.0/14.0;
static double SMALL_RF=5.0*DBL_MIN;
static double BIG_RF=0.2*DBL_MAX;

double RF(double x, double y, double z)
{
    if (MIN(MIN(x,y),z) < 0.0 ||
        MIN(MIN(x+y,x+z),y+z) < SMALL_RF ||
        MAX(MAX(x,y),z) > BIG_RF)
      error("RF(): Invalid arguments in function RF.");

    double A,sqrty,sqrtz,dx,dy,dz,lam;

    do {
        sqrty=sqrt(y);
        sqrtz=sqrt(z);
        lam=sqrt(x)*(sqrty+sqrtz)+sqrty*sqrtz;
        x=0.25*(x+lam);
        y=0.25*(y+lam);
        z=0.25*(z+lam);
        A=ONE_THIRD*(x+y+z);
        dx=(A-x)/A;
        dy=(A-y)/A;
        dz=(A-z)/A;
    } while (MAX(MAX(fabs(dx),fabs(dy)),fabs(dz)) > ERRTOL_RF);

    double e2=dx*dy-dz*dz;
    double e3=dx*dy*dz;

    return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(A);
}

double RD(double x, double y,double z)
{
    if (MIN(x,y) < 0.0 ||
        MIN(x+y,z) < SMALL_RD ||
        MAX(MAX(x,y),z) > BIG_RD)
      error("RD(): Invalid arguments in function RF.");


    double sum=0.0;
    double fac=1.0;
    double A,dx,dy,dz,sqrty,sqrtz,lam;

    do {
        sqrty=sqrt(y);
        sqrtz=sqrt(z);
        lam=sqrt(x)*(sqrty+sqrtz)+sqrty*sqrtz;
        sum += fac/(sqrtz*(z+lam));
        fac*=0.25;
        x=0.25*(x+lam);
        y=0.25*(y+lam);
        z=0.25*(z+lam);
        A=0.2*(x+y+3.0*z);
        dx=(A-x)/A;
        dy=(A-y)/A;
        dz=(A-z)/A;
    } while (MAX(MAX(fabs(dx),fabs(dy)),fabs(dz)) > ERRTOL_RD);

    double e1=dx*dy;
    double e2=dz*dz;
    double e3=e1-e2;
    double e4=e1-6.0*e2;
    double e5=e4+e3+e3;

    return 3.0*sum+fac*(1.0+e4*(-D1+D5*e4-D6*dz*e5)
           +dz*(D2*e5+dz*(-D3*e3+dz*D4*e1)))/(A*sqrt(A));

}

double elleptint(double phi,double k) {
  double sp=sin(phi);
  double qcp=SQR(cos(phi));
  double r=(1.0-sp*k)*(1.0+sp*k);

  return sp*(RF(qcp,r,1.0)-(SQR(sp*k))*RD(qcp,r,1.0)/3.0);
}
