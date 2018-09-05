/**
 * utils.cpp
 *
 *  Created on: 01.04.2014
 *      Author: franke
 */

#define BITMAP_GREY "P5"
#define BITMAP_RGB "P6"

#include "Utils.h"
#include <R_ext/Lapack.h>

/* show arguments of .External interface */
SEXP showArgs(SEXP args) {
  int i, nargs;
  Rcomplex cpl;
  const char *name;

  if((nargs = length(args) - 1) > 0) {
    for(i = 0; i < nargs; i++) {
      args = CDR(args);
      name = CHAR(PRINTNAME(TAG(args)));
      switch(TYPEOF(CAR(args))) {
      case REALSXP:
        Rprintf("[%d] '%s' %f\n", i+1, name, REAL(CAR(args))[0]);
        break;
      case LGLSXP:
      case INTSXP:
        Rprintf("[%d] '%s' %d\n", i+1, name, INTEGER(CAR(args))[0]);
        break;
      case CPLXSXP:
        cpl = COMPLEX(CAR(args))[0];
        Rprintf("[%d] '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
        break;
      case STRSXP:
        Rprintf("[%d] '%s' %s\n", i+1, name,
               CHAR(STRING_ELT(CAR(args), 0)));
        break;
      default:
        Rprintf("[%d] '%s' R type\n", i+1, name);
      }
    }
  }
  return(R_NilValue);
}

/* get a single call to an R function */
SEXP getCall(SEXP R_fname, SEXP R_args, SEXP R_rho) {
  SEXP RCallBack = R_NilValue;
  PROTECT(RCallBack = allocVector(LANGSXP, LENGTH(R_args)+1 ));
  SETCAR( RCallBack, findFun(install(CHAR(STRING_ELT(R_fname, 0))),R_rho ));

  SEXP p = CDR(RCallBack);
  SEXP names = R_NilValue;
  PROTECT(names = getAttrib(R_args, R_NamesSymbol));
  if(isNull(names))
   error(_("'R_NamesSymbol' is undefined in argument list."));

  /* debugging */
  //Rf_PrintValue(R_args);

  for (int i=0; p!=R_NilValue; p=CDR(p),i++) {
    SETCAR(p,VECTOR_ELT(R_args,i));
    SET_TAG(p,install(CHAR(STRING_ELT(names,i))));
  }

  UNPROTECT(2);
  return RCallBack;
}
/**
 * @brief Simple bivariate normal random variable
 */
void rbinorm(double mx, double sdx, double my, double sdy,double rho, double &x, double &y) {
  double q1=rnorm(0,1), q2=rnorm(0,1);
  x=std::sqrt(1-R_pow(rho,2.0))*sdx*q1 + rho*sdx*q2 + mx;
  y=my+sdy*q2;
}


/**
 * @brief Sample from x={0,1,2,3}
 *        according to given cumulative probabilities p
 *        Number of elements to sample from is n=4
 *
 * @param p         given cumulative probabilities
 */
int sample_k(double *p) {
  int j;
  double rU = unif_rand();
  for (j=0;j<3;j++) {
      if (rU <= p[j])
        break;
  }
  return j;
}

void rbinorm_exact(double *p, double mx, double sdx, double my, double sdy,
					 double rho, double &x, double &y)
{
  double q1=rnorm(0,1), q2=rnorm(0,1);

  int k = sample_k(p);
  x=std::sqrt(1-R_pow(rho,2.0))*sdx*q1 + rho*sdx*q2 + (mx+k*sdx*sdx);
  y=my+sdy*q2;
}

/**
 * @brief Calculate probabilities
 *
 * @param mx       mu of logN major semi-axis
 * @param sdx      sd of logN major semi-axis
 * @param lx       upper box limit x
 * @param ly       upper box limit y
 * @param lz       upper box limit z
 * @param p        [out] probabilities
 * @param mu       [out] factor for intensity
 */
void cum_prob_k(double mx, double sdx2, double lx, double ly, double lz, double *p, double *mu) {
  double a[4] = {lx*ly*lz,
                 2.0*(lx*ly+lx*lz+ly*lz),
                 M_PI*(lx+ly+lz),
                 4.0*M_PI/3.0};

  p[0] = a[0];
  p[1] = a[1]*std::exp(mx+0.5*sdx2);
  p[2] = a[2]*std::exp(2.0*(mx+sdx2));
  p[3] = a[3]*std::exp(3.0*mx+4.5*sdx2);

  double mk=0.0;
  for (int i=0;i<4; i++)
    mk += p[i];

  p[0] /= mk;
  p[1] /= mk;
  p[2] /= mk;
  p[3] /= mk;

  /* cumulative probabilities */
  for (int i=1;i<4; i++)
    p[i] += p[i-1];
  *mu=mk;
}


/* get the elements of a list */
SEXP getListElement(SEXP list, const char *str)
{
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;

}

SEXP getVar(SEXP name, SEXP rho)
{
    if(!isString(name) || length(name) != 1)
        error("name is not a single string");
    if(!isEnvironment(rho))
        error("rho should be an environment");

    return findVar(installChar(STRING_ELT(name, 0)), rho);
}
