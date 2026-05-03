#include "cgeneric_defs.h"
#include "stdio.h"

typedef struct 
{
    double *Q;
    double *lconst;
    double *theta;
}
my_cache_tp;

double *inla_cgeneric_rspde_1d_general_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {
  double *const_store = NULL;
  double*Q_store = NULL;
  double *ret = NULL;
  double ltau, lkappa, tau, kappa;
  char *prior_nu_dist;
  double nu;
  int N, M, i, k, rspde_order;
  double lnu;
  char *parameterization;

  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);
  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug == 1){
    debug = 1;
  }
  assert(!strcasecmp(data->doubles[0]->name, "nu_upper_bound"));
  double nu_upper_bound = data->doubles[0]->doubles[0];
  
  assert(!strcasecmp(data->mats[0]->name, "rational_table"));
  inla_cgeneric_mat_tp *rational_table = data->mats[0];
  assert(rational_table->nrow == 999);  
  
  assert(!strcasecmp(data->ints[2]->name, "graph_opt_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "graph_opt_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);
  
  assert(!strcasecmp(data->doubles[1]->name, "start_theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[1];

  assert(!strcasecmp(data->doubles[2]->name, "theta_prior_mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[2];
  
  assert(!strcasecmp(data->mats[1]->name, "theta_prior_prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[1];
  
  assert(!strcasecmp(data->doubles[3]->name, "prior_nu_loglocation"));
  double prior_nu_loglocation = data->doubles[3]->doubles[0];
  

  assert(!strcasecmp(data->doubles[4]->name, "prior_nu_mean"));
  double prior_nu_mean = data->doubles[4]->doubles[0];
 
  assert(!strcasecmp(data->doubles[5]->name, "prior_nu_prec"));
  double prior_nu_prec = data->doubles[5]->doubles[0];
  
  assert(!strcasecmp(data->doubles[6]->name, "prior_nu_logscale"));
  double prior_nu_logscale = data->doubles[6]->doubles[0];
  
  assert(!strcasecmp(data->doubles[7]->name, "start_nu"));
  double start_nu = data->doubles[7]->doubles[0];
 
  assert(!strcasecmp(data->chars[2]->name, "prior_nu_dist"));
  prior_nu_dist = &data->chars[2]->chars[0];
  
  
  assert(!strcasecmp(data->ints[4]->name, "rspde_order"));
  rspde_order = data->ints[4]->ints[0];

  assert(!strcasecmp(data->chars[3]->name, "parameterization"));
  parameterization = &data->chars[3]->chars[0];

  //assert(!strcasecmp(data->chars[4]->name, "prior.theta.param"));
  //char* theta_param = &data->chars[4]->chars[0];

  assert(!strcasecmp(data->doubles[8]->name, "loc"));
  double *loc = data->doubles[8]->doubles;
  int n_loc = data->doubles[8]->len;
  
  assert(!strcasecmp(data->ints[5]->name, "es"));
  int equally_spaced = data->ints[5]->ints[0];
  
  assert(!strcasecmp(data->ints[6]->name, "nu_fixed"));
  int nu_fixed = data->ints[6]->ints[0];

  int n_par = 3;
  
  
  if(nu_fixed == 1) {
      n_par = 2;
      nu = start_nu;
  } else {
      if (theta) {
          lnu = theta[n_par-1];
          nu = (exp(lnu)/(1.0 + exp(lnu))) * nu_upper_bound;
      }
      else {   
          lnu = nu = NAN;
      }    
  }
  
  
  //alpha = nu + 1.0 / 2.0;
  
  if (theta) {
    // interpretable parameters 
    if(!strcasecmp(parameterization, "matern")){
      lkappa = 0.5 * log(8.0 * nu) - theta[1];
      ltau = - theta[0] + 0.5 *(
        lgamma(nu) - 2.0 * nu * lkappa - 0.5 * log(4 * M_PI) - lgamma(nu + 0.5)
      );

    } else if(!strcasecmp(parameterization, "matern2")) {
      lkappa = - theta[1];
      ltau = - 0.5 * theta[0] + 0.5 *(
        lgamma(nu) - 2.0 * nu * lkappa - 0.5 * log(4 * M_PI) - lgamma(nu + 0.5)
      );
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
  
  double sigma = sqrt(tgamma(nu) / (pow(tau,2.0) * pow(kappa, 2.0*nu) * sqrt(4.0 * M_PI) * tgamma(nu + 0.5)));
  
  
  if (!(data->cache)) {
#pragma omp critical  (Name_7c3b4712ebb2dda8def3a5273e2a7e6cf1794b5d)
      if (!(data->cache)) {
          data->cache = (void **) Calloc(CGENERIC_CACHE_LEN(data), my_cache_tp *);
      }
  }
  
  int cache_idx;
  CGENERIC_CACHE_ASSIGN_IDX(cache_idx, data);
  my_cache_tp *cache = ((my_cache_tp **) data->cache)[cache_idx];
  
  
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
      
      double *r, *p, k_rat;
      if(rspde_order > 0) {
          int n_terms = 2*rspde_order + 2;
          double new_alpha = nu + 0.5;
          int row_nu = (int)round(1000*cut_decimals(new_alpha))-1;
          
          double *rat_coef = Calloc(n_terms-1, double);
          rat_coef = &rational_table->x[row_nu*n_terms+1];
          
          r = Calloc(rspde_order, double);
          p = Calloc(rspde_order, double);
          
          r = &rat_coef[0];
          p = &rat_coef[rspde_order];
          k_rat = rat_coef[2*rspde_order];    
      } else {
          r = Calloc(1, double);
          p = Calloc(1, double);
          
          r[0] = 1.0;
          p[0] = 0.0;
          k_rat = 1.0; 
      }
      
      const_store = Calloc(1, double);

      if(cache) {
          // cache exists, check if theta is the same
          
          if((n_par == 2 && cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) || 
             (n_par == 3 && cache->theta[0] == theta[0] && cache->theta[1] == theta[1] && cache->theta[2] == theta[2])) {

              memcpy(ret + 2, cache->Q, M * sizeof(double));

          } else {
              compute_Q1d(n_loc, loc, rspde_order, kappa, sigma, p, r, k_rat, &ret[k],
                          graph_i->ints, graph_j->ints, nu, M, equally_spaced, nu_upper_bound, N, const_store);
              
          }
      } else {
          compute_Q1d(n_loc, loc, rspde_order, kappa, sigma, p, r, k_rat, &ret[k],
                      graph_i->ints, graph_j->ints, nu, M, equally_spaced, nu_upper_bound, N, const_store);
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
      ret = Calloc(n_par+1, double);
      ret[0] = n_par;
      ret[1] = start_theta->doubles[0];
      ret[2] = start_theta->doubles[1];
      if(nu_fixed == 0) { 
          ret[n_par] = log(start_nu/(nu_upper_bound - start_nu));    
      } 
      
      break;
  }
      
  case INLA_CGENERIC_LOG_NORM_CONST:
  {
    
      ret = Calloc(1, double);
      Q_store = Calloc(M, double);
      
      double *r, *p, k_rat;
      if(rspde_order > 0) {
          int n_terms = 2*rspde_order + 2;
          double new_alpha = nu + 0.5;
          int row_nu = (int)round(1000*cut_decimals(new_alpha))-1;
          
          double *rat_coef = Calloc(n_terms-1, double);
          rat_coef = &rational_table->x[row_nu*n_terms+1];
          
          r = Calloc(rspde_order, double);
          p = Calloc(rspde_order, double);
          
          r = &rat_coef[0];
          p = &rat_coef[rspde_order];
          k_rat = rat_coef[2*rspde_order];    
      } else {
          r = Calloc(1, double);
          p = Calloc(1, double);
          
          r[0] = 1.0;
          p[0] = 0.0;
          k_rat = 1.0; 
      }
      
      
      if(cache) {
          // cache exists, check if theta is the same
          if((n_par == 2 && cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) || 
             (n_par == 3 && cache->theta[0] == theta[0] && cache->theta[1] == theta[1] && cache->theta[2] == theta[2])) {
              // same parameters, return const
              memcpy(ret, cache->lconst, sizeof(double));
          } else {
              // new parameters, compute the quantities
              compute_Q1d(n_loc, loc, rspde_order, kappa, sigma, p, r, k_rat, &Q_store[0],
                          graph_i->ints, graph_j->ints, nu, M, equally_spaced, nu_upper_bound, N, &ret[0]);  
              
              //then cache them and update theta
              memcpy(cache->Q, Q_store,  M * sizeof(double));
              memcpy(cache->lconst, ret, sizeof(double));
              memcpy(cache->theta, theta,  2 * sizeof(double));
          }
      } else {
          // cache empty, compute quantities
          compute_Q1d(n_loc, loc, rspde_order, kappa, sigma, p, r, k_rat, &Q_store[0],
                      graph_i->ints, graph_j->ints, nu, M, equally_spaced, nu_upper_bound, N, &ret[0]);  
          
          //then allocate memory in the cache and store it and theta
          ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
          cache->Q = Calloc(M, double);
          cache->lconst = Calloc(1, double);
          cache->theta = Calloc(2, double);
          memcpy(cache->Q, Q_store,  M * sizeof(double));
          memcpy(cache->lconst, ret,  sizeof(double));
          memcpy(cache->theta, theta,  2 * sizeof(double));
      }
      
      break;
  }
      
  case INLA_CGENERIC_LOG_PRIOR:
  {
      ret = Calloc(1, double);
      
      ret[0] = 0.0;
      
      if(nu_fixed == 0) {
          if(!strcasecmp(prior_nu_dist, "lognormal")){
              ret[0] += -0.5 * SQR(lnu - prior_nu_loglocation)/(SQR(prior_nu_logscale));
              ret[0] += -log(prior_nu_logscale) - 0.5 * log(2.0*M_PI);
              ret[0] -= log(pnorm(log(nu_upper_bound), prior_nu_loglocation, prior_nu_logscale));
              
          }
          else { 
              double s_1 = (prior_nu_mean / nu_upper_bound) * prior_nu_prec;
              double s_2 = (1 - prior_nu_mean / nu_upper_bound) * prior_nu_prec;
              ret[0] += logdbeta(nu / nu_upper_bound, s_1, s_2) - log(nu_upper_bound);
          }    
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