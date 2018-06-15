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

void *getExternalPtr(SEXP ext);
void checkPtr(SEXP ptr,SEXP type);
Rboolean isNullPtr(SEXP ptr, SEXP type);

SEXP getVar(SEXP name, SEXP rho);
SEXP getListElement (SEXP list, const char *str);

/* internally used */

#define GET_CALL(d,i) getCall( VECTOR_ELT(d->fname,i),VECTOR_ELT(d->args,i),d->rho)
#define GET_NAME(d,i) CHAR(STRING_ELT( VECTOR_ELT(d->fname,i), 0))


#ifdef __cplusplus
extern "C" {
#endif

/* pointer to some R's random generators */
typedef double (*rdist2_t)(double, double);

/* a const dummy function */
inline double rconst(double x, double dummy=0) { return x; }

/* bivariate normal random x,y*/
void rbinorm(double mx, double sdx, double my, double sdy,double rho, double &x, double &y);

/* exact bivariate normal random x,y : also samples integer k*/
void rbinorm_exact(double *p, double mx, double sdx, double my, double sdy,double rho, double &x, double &y);

/* sample */
void sample_k(double *p, int &k);

/* cumulative probabilities */
void cum_prob_k(double mx, double sdx2, double lx, double ly, double lz, double *p, double *mu);

/*
 * R call data struct, fname,args and call could also be
 * lists of names, args and calls to functions
 */
typedef struct R_Calldata_s {
    SEXP fname,args,rho,label,call;
    int nprotect, isPerfect;
} *R_Calldata;


#ifdef __cplusplus
}
#endif


template<typename R_TYPE>
struct R_eval_t {
  typedef R_TYPE return_type;
  SEXP call, rho;
  R_eval_t(SEXP _call, SEXP _rho) :  call(_call), rho(_rho)  {};
  inline return_type operator()() { return AS_NUMERIC(eval(call,rho));  }
};

template<>
struct R_eval_t<double> {
  SEXP call, rho;
  R_eval_t(SEXP _call, SEXP _rho) :  call(_call), rho(_rho)
  {
  };
  inline double operator()() {
    return asReal(eval(call,rho));
  }
};

struct R_rlnorm_t {
  double mx,sdx;
  R_rlnorm_t(double p,double q)
    : mx(p), sdx(q)
  {};

  inline double operator()() {  return rlnorm(mx,sdx); }
};

template<typename F >
struct R_rndGen_t {
  F rdist2;
  double mx,sdx;

  R_rndGen_t(double p,double q, const char* fname)
    : mx(p), sdx(q)
  {
    if ( !std::strcmp(fname, "rbeta" )) {
        rdist2=rbeta;
    } else if(!std::strcmp(fname, "rlnorm")) {
        rdist2=rlnorm;
    } else if(!std::strcmp(fname, "rgamma")) {
        rdist2=rgamma;
    } else if(!std::strcmp(fname, "runif" )) {
        rdist2=rweibull;
    } else if(!std::strcmp(fname, "const" )) {
        rdist2=rconst;
    } else {
        error("Undefined random generating function for radii distribution");
    }

  };

  inline double operator()() { return rdist2(mx,sdx); }
};




/**
 * \brief Show aruments of R function call
 *
 * @param args list of arguments
 * @return
 */
SEXP showArgs(SEXP args);

/**
 *
 * @param R_fname
 * @param R_arg
 * @param R_rho
 * @return
 */
SEXP getSingleCall(SEXP R_fname, SEXP R_arg, SEXP R_rho);

/**
 * \brief Define R call to user function from C level
 *
 * @param  R_fname (string) name of function to call
 * @param  R_args  list of arguments function
 * @param  R_rho   R environment
 * @return SEXP call as evaluated by eval(Call, Env)
 */
SEXP getCall(SEXP R_fname, SEXP R_args, SEXP R_rho);

/**
 * \brief Delete R_Calldata struct and and UNPROTECT SEXPs
 * @param call
 */
void deleteRCall(R_Calldata call);

R_Calldata getRCallParam(SEXP R_param, SEXP R_cond);

#endif /* UTILS_H_ */
