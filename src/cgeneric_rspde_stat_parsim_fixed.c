
#include "cgeneric_defs.h"
#include <stdio.h>


// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_stat_parsim_fixed_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {
  double *ret = NULL;
  double ltau, lkappa, tau, kappa, prior_theta1_meanlog, nu;
  double alpha;
  int m_alpha;
  double prior_theta1_sdlog, prior_theta2_meanlog, prior_theta2_sdlog;
  double start_theta1, start_theta2;
  int N, M, i, k;
  double d;
  char *parameterization;
  int fem_size;
  int one = 1;

   // the size of the model
  assert(data->n_ints == 4);

  // the number of doubles
  assert(data->n_doubles == 9);

  // the number of strings
  assert(data->n_chars == 3);

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

  assert(!strcasecmp(data->chars[2]->name, "parameterization"));
  parameterization = &data->chars[2]->chars[0];

  assert(!strcasecmp(data->doubles[0]->name, "d"));
  d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "nu"));
  nu = data->doubles[1]->doubles[0];

  alpha = nu + d / 2.0;
  m_alpha = (int) floor(alpha);

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

  assert(!strcasecmp(data->doubles[7]->name, "start.theta1"));
  start_theta1 = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "start.theta2"));
  start_theta2 = data->doubles[8]->doubles[0];

  if (theta) {
    // interpretable parameters 
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
    ltau = lkappa = tau = kappa = NAN;
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


        double *coeff;

        coeff = Calloc(m_alpha+2, double);

        coeff = markov_approx_coeff(alpha/2.0, kappa, (int)d);

      // FORTRAN IMPLEMENTATION

      dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
      coeff[0] = coeff[0] * SQR(tau);
      dscal_(&M, &coeff[0], &ret[k], &one);

      for(i = 1; i < m_alpha + 2; i++){
        coeff[i] = coeff[i] * SQR(tau);
        daxpy_(&M, &coeff[i], &fem->doubles[i*M], &one, &ret[k], &one);
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
      ret = Calloc(3, double);
      ret[0] = 2;
      ret[1] = start_theta1;
      ret[2] = start_theta2;
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
	  break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}