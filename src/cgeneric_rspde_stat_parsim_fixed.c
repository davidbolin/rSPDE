
#include "cgeneric_defs.h"
#include <stdio.h>


// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_stat_parsim_fixed_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {
  double *ret = NULL;
  double ltau, lkappa, tau, kappa, nu;
  double alpha;
  int m_alpha;
  int N, M, i, k;
  double d;
  char *parameterization, *theta_param;
  int fem_size;
  int one = 1;

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

  assert(!strcasecmp(data->chars[3]->name, "prior.theta.param"));
  theta_param = &data->chars[3]->chars[0];
  
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
  assert(!strcasecmp(data->doubles[3]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[3];

  assert(!strcasecmp(data->mats[0]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[0];

  assert(!strcasecmp(data->doubles[4]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[4];

  if (theta) {
    // interpretable parameters 
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
    }  else {
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
      ret[1] = start_theta->doubles[0];
      ret[2] = start_theta->doubles[1];
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

      if(!strcasecmp(theta_param, "theta") || !strcasecmp(parameterization, "spde")){
          ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                      theta_prior_prec->x, theta);
      }
      else {
        double theta_prior_mean_spde[2], theta_spde[2];
        theta_spde[1] = lkappa;
        theta_spde[0] = ltau;
        theta_prior_mean_spde[1] = 0.5 * log(8.0 * nu) - theta_prior_mean->doubles[1];
        theta_prior_mean_spde[0] = - theta_prior_mean->doubles[0] + 0.5 *(
          lgamma(nu) - 2.0 * nu * theta_prior_mean_spde[1] - 
          (d/2.0) * log(4 * M_PI) - lgamma(nu + d/2.0)
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