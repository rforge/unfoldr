/**
 * utils.h
 *
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <cstring>
#include <cmath>

#include "Rheaders.h"

/* functions with R types */

SEXP getVar(SEXP name, SEXP rho);
SEXP getListElement (SEXP list, const char *str);

/* internally used */

#define GET_CALL(fname,args,rho,i) getCall( VECTOR_ELT(fname,i),VECTOR_ELT(args,i),rho)
#define GET_NAME(fname,i) CHAR(STRING_ELT( VECTOR_ELT(fname,i), 0))


#ifdef __cplusplus
extern "C" {
#endif

/* pointer to some R's random generators */
typedef double (*rdist2_t)(double, double);

/* a const function */
inline double rconst(double x, double dummy=0) {  return (dummy=x); }

/* bivariate normal random x,y*/
void rbinorm(double mx, double sdx, double my, double sdy,double rho, double &x, double &y);

/* exact bivariate normal random x,y : also samples integer k*/
void rbinorm_exact(double *p, double mx, double sdx, double my, double sdy,double rho, double &x, double &y);

/* sample */
int sample_k(double *p);

/* cumulative probabilities */
void cum_prob_k(double mx, double sdx2, double lx, double ly, double lz, double *p, double *mu);

#ifdef __cplusplus
}
#endif


template<typename R_TYPE>
struct R_eval_t {
  SEXP call, rho;
  double mu;

  R_eval_t(SEXP _call, SEXP _rho, double _mu) :  call(_call), rho(_rho), mu(_mu)  {};

  R_TYPE operator()() { return eval(call,rho);  }
};


template<>
struct R_eval_t<double> {
  SEXP call, rho;
  double mu;

  R_eval_t(SEXP _call, SEXP _rho, double _mu) :  call(_call), rho(_rho), mu(_mu) {};

  double operator()() {
	  int info = 0;
	  SEXP Reval;

	  PROTECT(Reval = R_tryEval(call,rho,&info));
	  if(info != 0)
	    error(_("simJoint(): `try` error in user defined distribution function."));

	  UNPROTECT(1);
	  return REAL(Reval)[0];
  }
};


/**
 * \brief Show aruments of R function call
 *
 * @param args list of arguments
 * @return
 */
SEXP showArgs(SEXP args);

/**
 * \brief Define R call to user function from C level
 *
 * @param  R_fname (string) name of function to call
 * @param  R_args  list of arguments function
 * @param  R_rho   R environment
 * @return SEXP call as evaluated by eval(Call, Env)
 */
SEXP getCall(SEXP R_fname, SEXP R_args, SEXP R_rho);

#endif /* UTILS_H_ */
