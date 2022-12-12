#include "cgeneric_defs.h"
#include "stdio.h"
// #include "gsl/gsl_vector_double.h"


// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_nonstat_fixed_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

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


  assert(!strcasecmp(data->doubles[0]->name, "d"));
  double d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "r_ratapprox"));
  double *r = data->doubles[1]->doubles;

  assert(!strcasecmp(data->doubles[2]->name, "p_ratapprox"));
  double *p = data->doubles[2]->doubles;
  
  assert(!strcasecmp(data->doubles[3]->name, "k_ratapprox"));
  double k_rat = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "nu"));
  double nu = data->doubles[4]->doubles[0];

  double alpha = nu + d / 2.0;
  int m_alpha = floor(alpha);

  assert(!strcasecmp(data->smats[0]->name, "C"));
  inla_cgeneric_smat_tp *C = data->smats[0];
  
  assert(!strcasecmp(data->smats[1]->name, "G"));
  inla_cgeneric_smat_tp *G = data->smats[1];

  int n_mesh = C->ncol;

  assert(!strcasecmp(data->mats[0]->name, "B_tau"));
  inla_cgeneric_mat_tp *B_tau = data->mats[0];

  assert(!strcasecmp(data->mats[1]->name, "B_kappa"));
  inla_cgeneric_mat_tp *B_kappa = data->mats[1];

  int n_par = B_tau->ncol;

  // Starting values

  assert(!strcasecmp(data->doubles[5]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[5];
 
  assert(!strcasecmp(data->doubles[6]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[6];

  assert(!strcasecmp(data->mats[2]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[2];

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

      compute_Q_fixednu(n_mesh, C->x, C->i, C->j,
                        C->n,
                        G->x, G->i, G->j,
                        G->n,
                        B_kappa->x, B_tau->x,
                        B_kappa->ncol, rspde_order,
                        theta, p, r, k_rat,
                        m_alpha, &ret[k],
                        alpha);

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
      ret = Calloc(n_par, double);
      ret[0] = n_par-1;
      for(i=1; i<n_par; i++){
        ret[i] = start_theta->doubles[i-1];
      }
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