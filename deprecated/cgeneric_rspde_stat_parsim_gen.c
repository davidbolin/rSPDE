#include "cgeneric_defs.h"
#include <stdio.h>

double kappa_integral(int n, double beta, double kappa){
    double y;
    int k;
    y = 0;
    for(k = 0; k <= n; k++){
        y += (2*(k%2) - 1) * nChoosek(n,k)/(n-k-beta+1);
    }
    return(y*pow(kappa, 2*(n-beta+1)));
}


double *inla_cgeneric_rspde_stat_parsim_gen_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa, prior_theta1_meanlog, lnu, nu;
  double alpha, nu_upper_bound;
  int m_alpha;
  double prior_theta1_sdlog, prior_theta2_meanlog, prior_theta2_sdlog;
  double prior_nu_mean, prior_nu_loglocation, prior_nu_prec;
  double prior_nu_logscale;
  double start_theta1, start_theta2, start_nu;
  int N, M, i, k, j;
  double d;
  char *prior_nu_dist, *parameterization;
  int fem_size;
  int one = 1;
  double *Bmat, *coeff, *tmp_mat, *tmp_vec;
//   double *coeff;

   // the size of the model
  assert(data->n_ints == 4);

  // the number of doubles
  assert(data->n_doubles == 14);

  // the number of strings
  assert(data->n_chars == 4);

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

  assert(!strcasecmp(data->doubles[0]->name, "d"));
  d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "nu_upper_bound"));
  nu_upper_bound = data->doubles[1]->doubles[0];

  alpha = nu_upper_bound + d / 2.0;
  m_alpha = floor(alpha);

  Bmat = Calloc(SQR(m_alpha+2), double);
  coeff = Calloc(m_alpha+2, double);
  tmp_mat = Calloc(SQR(m_alpha+2), double);
  tmp_vec = Calloc(m_alpha+2, double);

  assert(!strcasecmp(data->doubles[2]->name, "matrices_full"));
  inla_cgeneric_vec_tp *fem = data->doubles[2];

  fem_size = (fem->len)/(m_alpha+2);
  assert(M == fem_size);

  // prior parameters
  assert(!strcasecmp(data->doubles[3]->name, "prior.theta1.meanlog"));
  prior_theta1_meanlog = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "prior.theta1.sdlog"));
  prior_theta1_sdlog = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior.theta2.meanlog"));
  prior_theta2_meanlog = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "prior.theta2.sdlog"));
  prior_theta2_sdlog = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "prior.nu.loglocation"));
  prior_nu_loglocation = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "prior.nu.mean"));
  prior_nu_mean = data->doubles[8]->doubles[0];

  assert(!strcasecmp(data->doubles[9]->name, "prior.nu.prec"));
  prior_nu_prec = data->doubles[9]->doubles[0];

  assert(!strcasecmp(data->doubles[10]->name, "prior.nu.logscale"));
  prior_nu_logscale = data->doubles[10]->doubles[0];

  assert(!strcasecmp(data->doubles[11]->name, "start.theta1"));
  start_theta1 = data->doubles[11]->doubles[0];

  assert(!strcasecmp(data->doubles[12]->name, "start.theta2"));
  start_theta2 = data->doubles[12]->doubles[0];

  assert(!strcasecmp(data->doubles[13]->name, "start.nu"));
  start_nu = data->doubles[13]->doubles[0];

  if (theta) {
    // interpretable parameters 
    lnu = theta[2];
    nu = (exp(lnu)/(1.0 + exp(lnu))) * nu_upper_bound;
    if(!strcasecmp(parameterization, "matern")){
      ltau = - theta[0];
      lkappa = 0.5 * log(8.0 * nu) - theta[1];
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

    //   coeff = (double *)malloc((new_m_alpha+2)*sizeof(double));

    int size = new_m_alpha + 2;
    double L = new_alpha - floor(new_alpha);
    for (i = 0; i < size; i++){
        coeff[i] = 2*kappa_integral(i,-new_alpha+2.0*(size-1)+1.0+L,kappa);
    for (j = 0; j < size; j++){
      Bmat[j+i*size] = 2*kappa_integral(i+j,2.0*(size-1)+1.0+L,kappa);
    }
  }


    int ipivot[size];
    int info;

    int sizep1 = size+1;

    dcopy_(&size, Bmat, &sizep1, tmp_vec, &one);

    for(i = 0; i < size; i++){
        tmp_vec[i] = sqrt(tmp_vec[i]);
    }


    dcopy_(&size, tmp_vec, &one, tmp_mat, &sizep1);

    dgesv_(&size, &one, tmp_mat, &size, ipivot, coeff, &size, &info);
    
    dgesv_(&size, &size, tmp_mat, &size, ipivot, Bmat, &size, &info);

    dgesv_(&size, &one, Bmat, &size, ipivot, coeff, &size, &info);

  assert(info == 0);

  double fact = exp(lgamma(nu))/(signgam*exp(lgamma(new_alpha))*pow((4.0*M_PI),d/2.0)*pow(kappa,(2*nu)));
    for(i = 0; i < size; i++){
        coeff[i] *= SQR(tau)*fact;
    }

      printf("START TO PRINT!\n");
      for(i = 0; i< new_m_alpha+2; i++){
        printf("coeff[%d] = %f\n", i, coeff[i]/SQR(tau));
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
      ret[1] = start_theta1;
      ret[2] = start_theta2;
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

      ret[0] += -0.5 * SQR(theta[0] - prior_theta1_meanlog)/(SQR(prior_theta1_sdlog)) - 
      log(prior_theta1_sdlog) - 0.5 * log(2.0 * M_PI);

      ret[0] += -0.5 * SQR(theta[1] - prior_theta2_meanlog)/(SQR(prior_theta2_sdlog)) - 
      log(prior_theta2_sdlog) - 0.5 * log(2.0 * M_PI);

      if(!strcasecmp(prior_nu_dist, "lognormal")){
        ret[0] += -log(nu) -0.5 * SQR(lnu - prior_nu_loglocation)/(SQR(prior_nu_logscale));
        ret[0] += -log(prior_nu_logscale) - 0.5 * log(2.0*M_PI);
        ret[0] -= log(pnorm(log(nu_upper_bound), prior_nu_loglocation, prior_nu_logscale));
      }
      else { // if(!strcasecmp(prior_nu_dist, "beta")){
        double s_1 = (prior_nu_mean / nu_upper_bound) * prior_nu_prec;
        double s_2 = (1 - prior_nu_mean / nu_upper_bound) * prior_nu_prec;
        ret[0] += logdbeta(nu / nu_upper_bound, s_1, s_2) - log(nu_upper_bound);
      }

	  break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}