#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))


# ifdef __SUPPORT_SNAN__
#  define iszero(x) (fpclassify (x) == FP_ZERO)
# else
#  define iszero(x) (((__typeof (x)) (x)) == 0)
# endif

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif
#define F_ONE ((fortran_charlen_t)1)


// https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c

unsigned nChoosek( int n, int k );
double cut_decimals(double nu);

void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);

void dscal_(int* N, double* DA, double* DX,int* INCX);

void dcopy_(int* N, double* DX, int* INCX, double* DY,int* INCY);

void daxpby_(int* N, double* DA, double* DX, int* INCX, double* DB, double* DY, int* INCY, double* DZ);

inla_cgeneric_func_tp inla_cgeneric_rspde_stat_int_model;

