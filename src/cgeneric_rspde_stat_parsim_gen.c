#include "cgeneric_defs.h"
// #include <stdio.h>

// Creates a diagonal matrix with the elements of the vector as diagonal entries
void diag(double *vec_long, double *vec_short, int n){ //n = size of vec_short
    int nplusone = n+1, one = 1;
    dcopy_(&n, vec_short, &one, vec_long, &nplusone);
}

// Mat should be given by column!
int solveAb(double *mat, double *vec, int size){  // compute A^{-1}b, where b is a vector
    int one=1;
    int ipivot[size];
    int info;
    dgesv_(&size, &one, mat, &size, ipivot, vec, &size, &info);
    return(info);
}

//Get diagonal of a square matrix
// n is the size (number of columns) of the matrix
void getDiag(double *mat, double *destvec, int n){
    int nplusone = n+1, one = 1;
    dcopy_(&n, mat, &nplusone, destvec, &one);
}

// Computes solve(solve(diag(sqrt(diag(B))),B),solve(diag(sqrt(diag(B))),c))
// This will also change the matrix mat (but for our application there is no problem)!
int CrazySolve(double *mat, double *in_out_vec, int size){
    double *tmp_vec, *tmp_mat;
    int i, ipivot[size];
    tmp_vec = Calloc(size, double);
    int info;
    getDiag(mat, tmp_vec, size);
    for(i = 0; i < size; i++){
        tmp_vec[i] = sqrt(tmp_vec[i]);
    }
    tmp_mat = Calloc(size*size, double);
    diag(tmp_mat, tmp_vec, size);
    solveAb(tmp_mat, in_out_vec, size);
    
    dgesv_(&size, &size, tmp_mat, &size, ipivot, mat, &size, &info);

    solveAb(mat, in_out_vec, size);

    free(tmp_vec);
    free(tmp_mat);
    return(info);
} 

double kappa_integral(int n, double beta, double kappa){
    double y;
    int k;
    y = 0;
    for(k = 0; k <= n; k++){
        y += (2*(k%2) - 1) * nChoosek(n,k)/(n-k-beta+1);
    }
    return(y*pow(kappa, 2*(n-beta+1)));
}

double * markov_approx_coeff(double beta, double kappa, int d){
    double nu = 2*beta - d/2.0;
    double alpha = nu + d/2.0;
    double L = alpha - floor(alpha);
    int p = (int)ceil(alpha);
    int i,j;
    double *Bmat;
    Bmat = Calloc( SQR(p+1), double);
    double *c_vec;
    c_vec = Calloc(p+1, double);
    for(i = 0; i <= p; i++){
        c_vec[i] = 1.0;
    }
    for (i = 0; i <= p; i++){
    c_vec[i] = 2*kappa_integral(i,-alpha+2.0*p+1.0+L,kappa);
    for (j = 0; j <= p; j++){
      Bmat[j+i*(p+1)] = 2*kappa_integral(i+j,2.0*p+1.0+L,kappa);
    }
  }
  int info;
  info = CrazySolve(Bmat, c_vec, p+1);
  assert(info == 0);
  // double fact = exp(lgamma(nu))/(signgam*exp(lgamma(alpha))*pow((4.0*M_PI),d/2.0)*pow(kappa,(2*nu)));
  //   for(i = 0; i <= p; i++){
  //       c_vec[i] *= fact;
  //   }
    return(c_vec);
}


double *inla_cgeneric_rspde_stat_parsim_gen_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa, lnu, nu;
  double alpha, nu_upper_bound;
  int m_alpha;
  double prior_nu_mean, prior_nu_loglocation, prior_nu_prec;
  double prior_nu_logscale;
  double start_nu;
  int N, M, i, k, j;
  double d;
  char *prior_nu_dist, *parameterization, *theta_param;
  int fem_size;
  int one = 1;
//   double *coeff;


  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug == 1){
    debug = 1;
  }

  assert(!strcasecmp(data->ints[2]->name, "graph_opt_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "graph_opt_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);

  assert(!strcasecmp(data->chars[2]->name, "prior.nu.dist"));
  prior_nu_dist = &data->chars[2]->chars[0];

  assert(!strcasecmp(data->chars[3]->name, "parameterization"));
  parameterization = &data->chars[3]->chars[0];

  assert(!strcasecmp(data->chars[4]->name, "prior.theta.param"));
  theta_param = &data->chars[4]->chars[0];

  assert(!strcasecmp(data->doubles[0]->name, "d"));
  d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "nu.upper.bound"));
  nu_upper_bound = data->doubles[1]->doubles[0];

  alpha = nu_upper_bound + d / 2.0;
  m_alpha = floor(alpha);

  assert(!strcasecmp(data->doubles[2]->name, "matrices_full"));
  inla_cgeneric_vec_tp *fem = data->doubles[2];

  fem_size = (fem->len)/(m_alpha+2);
  assert(M == fem_size);

  // prior parameters

  assert(!strcasecmp(data->doubles[3]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[3];

  assert(!strcasecmp(data->mats[0]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[0];

  assert(!strcasecmp(data->doubles[4]->name, "prior.nu.loglocation"));
  prior_nu_loglocation = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior.nu.mean"));
  prior_nu_mean = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "prior.nu.prec"));
  prior_nu_prec = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "prior.nu.logscale"));
  prior_nu_logscale = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[8];

  assert(!strcasecmp(data->doubles[9]->name, "start.nu"));
  start_nu = data->doubles[9]->doubles[0];

  if (theta) {
    // interpretable parameters 
    lnu = theta[2];
    nu = (exp(lnu)/(1.0 + exp(lnu))) * nu_upper_bound;
    if(!strcasecmp(parameterization, "matern")){
      lkappa = 0.5 * log(8.0 * nu) - theta[1];
      ltau = - theta[0] + 0.5 *(
        lgamma(nu) - 2.0 * nu * lkappa - (d/2.0) * log(4 * M_PI) - lgamma(nu + d/2.0)
      );
    } else if(!strcasecmp(parameterization, "matern2")) {
      lkappa = - theta[1];
      ltau = - 0.5 * theta[0] + 0.5 *(
        lgamma(nu) - 2.0 * nu * lkappa - (d/2.0) * log(4 * M_PI) - lgamma(nu + d/2.0)
      );
    } else {
      ltau = theta[0];
      lkappa = theta[1];
    }
    tau = exp(ltau);
    kappa = exp(lkappa);
  }
  else {   
    ltau = lkappa = lnu = tau = kappa = nu = NAN;
  }
  
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }
    
  case INLA_CGENERIC_GRAPH:
    {
      k=2;
      ret = Calloc(k + 2 * M, double);
      ret[0] = N;        	       /* dimension */
      ret[1] = M;		   /* number of (i <= j) */
      for (i = 0; i < M; i++) {
	      ret[k++] = graph_i->ints[i];
      }
      for (i = 0; i < M; i++) {
	      ret[k++] = graph_j->ints[i];
      }
      break;
    }
    
  case INLA_CGENERIC_Q:
    {
      k = 2;
      ret = Calloc(k + M, double);
      ret[0] = -1;		/* REQUIRED */
      ret[1] = M;		/* REQUIRED */

      double new_alpha = nu + d / 2.0;
      int new_m_alpha = (int) floor(new_alpha);

      if(new_alpha / 2.0 == (int) new_alpha/2.0){
      double sqkappatau1 = SQR(tau) * pow(kappa, 2 * new_m_alpha);
      double sqkappatau2 = SQR(tau) * new_m_alpha * pow(kappa, 2 * (new_m_alpha - 1));
      double sqtaukappatmp;
      dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
        dscal_(&M, &sqkappatau1, &ret[k], &one);
        daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
        if(new_m_alpha>=2){
          for(j = 2; j<= new_m_alpha; j++){
            sqtaukappatmp = SQR(tau) * pow(kappa, 2*(new_m_alpha-j)) * nChoosek(new_m_alpha, j);
            daxpy_(&M, &sqtaukappatmp, &fem->doubles[j*M], &one, &ret[k], &one);
          }
        }
      } else {

        double *coeff;

        coeff = Calloc(new_m_alpha+2, double);

        coeff = markov_approx_coeff(new_alpha/2.0, kappa, (int)d);
        
        for(i = 0; i < new_m_alpha + 2; i++){
            coeff[i] *= SQR(tau);
        }

      // FORTRAN IMPLEMENTATION

      dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
      dscal_(&M, &coeff[0], &ret[k], &one);

      for(i = 1; i < new_m_alpha + 2; i++){
        daxpy_(&M, &coeff[i], &fem->doubles[i*M], &one, &ret[k], &one);
      }
    //   free(coeff);
      }


      break;
    }
    
  case INLA_CGENERIC_MU:
    {
      ret = Calloc(1, double);
      ret[0] = 0.0;
      break;
    }
    
  case INLA_CGENERIC_INITIAL:
    {
      // return c(P, initials)
      // where P is the number of hyperparameters      
      ret = Calloc(4, double);
      ret[0] = 3;
      ret[1] = start_theta->doubles[0];
      ret[2] = start_theta->doubles[1];
      ret[3] = log(start_nu/(nu_upper_bound - start_nu));
      break;
    }
    
  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      break;
    }
    
  case INLA_CGENERIC_LOG_PRIOR:
    {
      ret = Calloc(1, double);

      ret[0] = 0.0;

      if(!strcasecmp(prior_nu_dist, "lognormal")){
        ret[0] += -0.5 * SQR(lnu - prior_nu_loglocation)/(SQR(prior_nu_logscale));
        ret[0] += -log(prior_nu_logscale) - 0.5 * log(2.0*M_PI);
        ret[0] -= log(pnorm(log(nu_upper_bound), prior_nu_loglocation, prior_nu_logscale));
      }
      else { // if(!strcasecmp(prior_nu_dist, "beta")){
        double s_1 = (prior_nu_mean / nu_upper_bound) * prior_nu_prec;
        double s_2 = (1 - prior_nu_mean / nu_upper_bound) * prior_nu_prec;
        ret[0] += logdbeta(nu / nu_upper_bound, s_1, s_2) - log(nu_upper_bound);
      }

      if(!strcasecmp(theta_param, "theta") || !strcasecmp(parameterization, "spde")){
          ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                      theta_prior_prec->x, theta);
      }
      else {
        double theta_prior_mean_spde[2], theta_spde[2], prior_nu_tmp;
        if(!strcasecmp(prior_nu_dist, "lognormal")){
          prior_nu_tmp = exp(prior_nu_loglocation);
        }
        else{
          prior_nu_tmp = prior_nu_mean;
        }
        theta_spde[1] = lkappa;
        theta_spde[0] = ltau;
        theta_prior_mean_spde[1] = 0.5 * log(8.0 * prior_nu_tmp) - theta_prior_mean->doubles[1];
        theta_prior_mean_spde[0] = - theta_prior_mean->doubles[0] + 0.5 *(
          lgamma(prior_nu_tmp) - 2.0 * prior_nu_tmp * theta_prior_mean_spde[1] - 
          (d/2.0) * log(4 * M_PI) - lgamma(prior_nu_tmp + d/2.0)
        );

        ret[0] += logmultnormvdens(2, theta_prior_mean_spde,
                                      theta_prior_prec->x, theta_spde);
      }

	  break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}