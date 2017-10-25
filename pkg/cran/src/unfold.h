/**
 * discretization.h
 *
 *  Created on: 18.03.2015
 *      Author: franke
 */

#ifndef BINNING_H_
#define BINNING_H_

#include "Rheaders.h"

template<class T>
inline T cot(const T x) {return tan(M_PI_2 - x);}

#ifdef __cplusplus
extern "C" {
#endif

  /* DLL Init */
  void R_init_unfoldr(DllInfo *info);

  /* .Call methods  */
  SEXP Binning1d(SEXP Rx, SEXP Rbin_x);
  SEXP Binning3d(SEXP Rx, SEXP Ry, SEXP Rz, SEXP Rbin_x, SEXP Rbin_y, SEXP Rbin_z);
  SEXP EMS(SEXP R_P, SEXP R_F, SEXP R_cond);
  SEXP CoefficientMatrixSpheroids(SEXP R_C, SEXP R_alpha, SEXP R_S, SEXP R_c, SEXP R_Theta, SEXP R_s, SEXP R_args );

  /* .C methods */
  void extern_findIndex(double *x, double *v, int *n, int *idx);
  void em_saltykov(int *_n, int *_nla, double *p, double *y, double *theta);
  void em_saltykov_p(int *nn, double *bin, double *xp);

#ifdef __cplusplus
}
#endif


#endif /* BINNING_H_ */
