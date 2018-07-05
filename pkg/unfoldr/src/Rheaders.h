/**
 * Rheaders.h
 *
 */

#ifndef RHEADERS_H_
#define RHEADERS_H_

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

/* definitions not involving SEXPs, plus _() */

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("unfoldr", String)
#else
#define _(String) (String)
#endif

/////////////////// Some R convenience macros ////////////////////////////////////

#define getDims(A) INTEGER(coerceVector (getAttrib(A, R_DimSymbol ) , INTSXP) )
#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define REAL_ARG_LIST(A,i) REAL(VECTOR_ELT( (A), (i) ))[0]


#define SET_CLASS_NAME(RObject,ClassName) {            \
  SEXP RClass;                                         \
  PROTECT(RClass=allocVector(STRSXP,1));               \
  SET_STRING_ELT( (RClass) ,0,mkChar( (ClassName) ));  \
  classgets( (RObject), (RClass) );                    \
  UNPROTECT(1);                                        \
}

#endif /* RHEADERS_H_ */
