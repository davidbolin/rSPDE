#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

// https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c

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
                    double *Q_out,
                    int *graph_i, int *graph_j, int M,
                    int matern_par, double start_nu, double nu, double d);

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


void compute_Q_alpha2(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int *lower_edges,
                                        int *upper_edges, int lower_edges_len, int upper_edges_len);

void compute_Q_alpha1_directional(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int stat_ind_len, int *stat_ind, int BC);


double adjusted_inv_logit(double z, double L);
double logit_adjusted(double x, double L);

double forward_nu(double lnu, double nu_upper_bound);

// Inverse transformation: Compute lnu from nu
double inverse_lnu(double nu, double nu_upper_bound);

#ifdef __cplusplus
extern "C" {
#endif

double cut_decimals(double nu);

double nChoosek(int n, int k);
void compute_Q_intrinsic(int size, 
                         double *entries_C, int *i_C, int *j_C, int n_nonzero_C,
                         double *entries_G, int *i_G, int *j_G, int n_nonzero_G,
                         double *theta, double *Q_out, int alpha,
                         int compute_Q, int compute_mean, int compute_logdet,
                         double *ld_out, double *mean_out);

void compute_Q_dim1(
    double kappa, double sigma, double gamma, double rho, int beta, int alpha,
    double* result,
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list);

void compute_Q_dim2(
    double kappa, double sigma, double gamma, double rho_1, double rho_2, int beta, int alpha,
    double* result,
    inla_cgeneric_smat_tp** Gtlist, 
    inla_cgeneric_smat_tp** Ctlist, 
    inla_cgeneric_smat_tp** B0list, 
    inla_cgeneric_smat_tp*** M2list);

void compute_Q_anisotropic(
    double hx, double hy, double hxy, double sigma, double nu,
    const inla_cgeneric_smat_tp *C,
    const inla_cgeneric_smat_tp *Ci,
    const inla_cgeneric_smat_tp *Hxx,
    const inla_cgeneric_smat_tp *Hyy,
    const inla_cgeneric_smat_tp *Hxy,
    const inla_cgeneric_smat_tp *Q_graph,    
    double *ret, int rspde_order,
    const inla_cgeneric_mat_tp *rational_table,
    int est_nu);

void compute_Q1d(int n, double *loc, int rspde_order, double kappa,
                 double sigma, double *rat_p, double *rat_r, double rat_k,
                 double *Q_out, int *graph_i, int *graph_j, double nu, int M,
                 int equally_spaced, double nu_upper_bound, int N, double *lconst);

void compute_Q_fintrinsic(double tau, double nu, 
                          const inla_cgeneric_smat_tp *C,
                          const inla_cgeneric_smat_tp *Ci,
                          const inla_cgeneric_smat_tp *G,
                          const inla_cgeneric_smat_tp *Q_graph,    
                          double *ret, int rspde_order,
                          const inla_cgeneric_mat_tp *rational_table,
                          int est_nu, int d, int compute_Q, 
                          int compute_mean, int compute_logdet,
                          double*ld_out, double *mean_out,
                          double scaling,
                          const inla_cgeneric_smat_tp *D);

/* ----- hybrid Whittle-Matern (alpha = 2) helpers ------------------ */

/* Fill Q_out with the M lower-triangular entries of
 *   Q = tau^2 * (G + kappa^2 C) * C^{-1} * (G + kappa^2 C)
 * at the (graph_i, graph_j) positions, where C is the mass-lumped FEM
 * mass matrix (passed in via its diagonal C_diag of length N). G is a
 * sparse triplet (rows G_i, cols G_j, values G_x, length G_n).
 */
void compute_Q_hybrid_alpha2(int N,
                             const double *C_diag,
                             const int *G_i, const int *G_j,
                             const double *G_x, int G_n,
                             double tau, double kappa,
                             const int *graph_i, const int *graph_j, int M,
                             double *Q_out);

/* Fill mu_out (length N) with mu = (G + kappa_mu^2 C)^{-1} (C * X * beta_X)
 * where X is a dense N x p matrix stored row by row and beta_X has
 * length p. C_diag is the mass-lumped diagonal of length N. The
 * `kappa_mu` argument lets the operator applied to the mean use a
 * range parameter that may differ from the kappa of the covariance.
 */
void compute_mu_hybrid_alpha2(int N,
                              const double *C_diag,
                              const int *G_i, const int *G_j,
                              const double *G_x, int G_n,
                              double kappa_mu,
                              const double *X_x, int p,
                              const double *beta_x,
                              double *mu_out);

#ifdef __cplusplus
}
#endif