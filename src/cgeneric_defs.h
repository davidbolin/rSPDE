#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

// https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c

double nChoosek( int n, int k );
double cut_decimals(double nu);

void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);

void dscal_(int* N, double* DA, double* DX,int* INCX);

void dcopy_(int* N, double* DX, int* INCX, double* DY,int* INCY);

void daxpby_(int* N, double* DA, double* DX, int* INCX, double* DB, double* DY, int* INCY, double* DZ);

void dgesv_(int *n, int *nrhs,  double *a,  int  *lda,  
           int *ipivot, double *b, int *ldb, int *info) ;

void dgemv_(char* trans, int* M, int* N, double* alpha, double* A,
           int* LDA, double* x, int* incx,
           double* beta, double* y, int* inc);

double * markov_approx_coeff(double beta, double kappa, int d);

double pnorm(double x, double mu, double sd);

double logdbeta(double x, double s_1, double s_2);

double logmultnormvdens(int npar, double *entries_mean, 
                        double *entries_prec,
                        double *entries_val);

void compute_Q(int size, double *entries_C, int *i_C, int *j_C,
                    int n_nonzero_C,
                    double *entries_G, int *i_G, int *j_G,
                    int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, int rspde_order, double *theta_entries,
                    double *rat_p, double *rat_r, double rat_k,
                    int m_alpha, double *Q_out,
                    int *graph_i, int *graph_j, int M, double alpha);

void compute_Q_fixednu(int size, double *entries_C, int *i_C, int *j_C,
                        int n_nonzero_C,
                        double *entries_G, int *i_G, int *j_G,
                        int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, int rspde_order, double *theta_entries,
                    double *rat_p, double *rat_r, double rat_k,
                    int m_alpha, double *Q_out, double alpha);

void compute_Q_integer(int size, double *entries_C, int *i_C, int *j_C,
                        int n_nonzero_C,
                        double *entries_G, int *i_G, int *j_G,
                        int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, double *theta_entries,
                    double *Q_out, int alpha);