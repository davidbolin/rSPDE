#include "cgeneric_defs.h"
#include "stdio.h"
// #include "gsl/gsl_vector_double.h"


// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_nonstat_general_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  char *prior_nu_dist;

  int k, i;

  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  int N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug == 1){
    debug = 1;
  }

  assert(!strcasecmp(data->ints[2]->name, "graph_opt_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  int M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "graph_opt_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);

  assert(!strcasecmp(data->ints[4]->name, "rspde_order"));
  int rspde_order = data->ints[4]->ints[0];

  assert(!strcasecmp(data->ints[5]->name, "matern_par"));
  int matern_par = data->ints[5]->ints[0];


  assert(!strcasecmp(data->doubles[0]->name, "d"));
  double d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "nu_upper_bound"));
  double nu_upper_bound = data->doubles[1]->doubles[0];

  assert(!strcasecmp(data->mats[0]->name, "rational_table"));
  inla_cgeneric_mat_tp *rational_table = data->mats[0];
  assert(rational_table->nrow == 999);  

  assert(!strcasecmp(data->smats[0]->name, "C"));
  inla_cgeneric_smat_tp *C = data->smats[0];
  
  assert(!strcasecmp(data->smats[1]->name, "G"));
  inla_cgeneric_smat_tp *G = data->smats[1];

  int n_mesh = C->ncol;

  assert(!strcasecmp(data->mats[1]->name, "B_tau"));
  inla_cgeneric_mat_tp *B_tau = data->mats[1];

  assert(!strcasecmp(data->mats[2]->name, "B_kappa"));
  inla_cgeneric_mat_tp *B_kappa = data->mats[2];

  int n_par = B_tau->ncol;

  // Prior param

  assert(!strcasecmp(data->doubles[2]->name, "prior.nu.loglocation"));
  double prior_nu_loglocation = data->doubles[2]->doubles[0];

  assert(!strcasecmp(data->doubles[3]->name, "prior.nu.logscale"));
  double prior_nu_logscale = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "prior.nu.mean"));
  double prior_nu_mean = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior.nu.prec"));
  double prior_nu_prec = data->doubles[5]->doubles[0];

  // Nu prior

  assert(!strcasecmp(data->chars[2]->name, "prior.nu.dist"));
  prior_nu_dist = &data->chars[2]->chars[0];

  // Starting values

  assert(!strcasecmp(data->doubles[6]->name, "start.nu"));
  double start_nu = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[7];
  
  assert(!strcasecmp(data->doubles[8]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[8];

  assert(!strcasecmp(data->mats[3]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[3];
  
  double lnu, nu;

   if (theta) {
    // interpretable parameters 
    lnu = theta[n_par-1];
    nu = (exp(lnu)/(1.0 + exp(lnu))) * nu_upper_bound;
  }
  else {   
    lnu = nu = NAN;
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

      int n_terms = 2*rspde_order + 2;

      double new_alpha = nu + d/2.0;

      int row_nu = (int)round(1000*cut_decimals(new_alpha))-1;

      double *rat_coef = Calloc(n_terms-1, double);
      
      rat_coef = &rational_table->x[row_nu*n_terms+1];

      double *r, *p, k_rat;

      r = Calloc(rspde_order, double);
      p = Calloc(rspde_order, double);
      
      r = &rat_coef[0];
      p = &rat_coef[rspde_order];
      k_rat = rat_coef[2*rspde_order];

      compute_Q(n_mesh, C->x, C->i, C->j,
                        C->n,
                        G->x, G->i, G->j,
                        G->n,
                        B_kappa->x, B_tau->x,
                        B_kappa->ncol, rspde_order,
                        theta, p, r, k_rat, &ret[k],
                        graph_i->ints, graph_j->ints,
                        M, matern_par, start_nu, nu,
                        d);

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
      ret = Calloc(n_par+1, double);
      ret[0] = n_par;
      for(i=1; i<n_par; i++){
        ret[i] = start_theta->doubles[i-1];
      }
      ret[n_par] = log(start_nu/(nu_upper_bound - start_nu));
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

      ret[0] += logmultnormvdens(n_par-1, theta_prior_mean->doubles,
                                  theta_prior_prec->x, theta);


	  break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}