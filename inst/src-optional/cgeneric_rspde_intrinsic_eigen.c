#include "cgeneric_defs.h"
#include "cgeneric.h"
#include "stdio.h"
// #include "gsl/gsl_vector_double.h"


// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_intrinsic_eigen(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double *mu_store = NULL;
  double *const_store = NULL;
  double*Q_store = NULL;
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

  assert(!strcasecmp(data->smats[0]->name, "C"));
  inla_cgeneric_smat_tp *C = data->smats[0];
  
  assert(!strcasecmp(data->smats[1]->name, "G"));
  inla_cgeneric_smat_tp *G = data->smats[1];
  
  assert(!strcasecmp(data->doubles[0]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[0];
  
  assert(!strcasecmp(data->mats[0]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[0];
  
  assert(!strcasecmp(data->doubles[1]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[1];
  
  assert(!strcasecmp(data->ints[4]->name, "alpha"));
  int alpha = data->ints[4]->ints[0];

  assert(!strcasecmp(data->ints[5]->name, "mean_correction"));
  int mean_correction = data->ints[5]->ints[0];
  
  assert(!strcasecmp(data->ints[6]->name, "true_scaling"));
  int true_scaling = data->ints[6]->ints[0];
  
  int n_mesh = C->ncol;
  
  
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

      
      compute_Q_intrinsic(n_mesh, 
                          C->x, C->i, C->j, C->n,
                          G->x, G->i, G->j, G->n,
                          theta, 
                          &ret[k],
                          alpha,
                          1,
                          0,
                          0,
                          const_store, 
                          mu_store);    

      break;
    }
    
  case INLA_CGENERIC_MU:
    {
      if(mean_correction == 1) { 
        ret = Calloc(1 + n_mesh, double);
        ret[0] = n_mesh;		/* REQUIRED */
        
        compute_Q_intrinsic(n_mesh, 
                            C->x, C->i, C->j, C->n,
                            G->x, G->i, G->j, G->n,
                            theta, 
                            Q_store,
                            alpha,
                            0,
                            1,
                            0,
                            const_store, 
                            &ret[1]);    
      } else {
          ret = Calloc(1, double);
          ret[0] = 0.0;
      }
      break;
    }
    
  case INLA_CGENERIC_INITIAL:
    {
      // return c(P, initials)
      // where P is the number of hyperparameters      
      ret = Calloc(2, double);
      ret[0] = 2;
      ret[1] = start_theta->doubles[0];
      ret[2] = start_theta->doubles[1];
      break;
    }
    
  case INLA_CGENERIC_LOG_NORM_CONST:
    {
        if(true_scaling) {
            ret = Calloc(1, double);
            compute_Q_intrinsic(n_mesh, 
                                  C->x, C->i, C->j, C->n,
                                  G->x, G->i, G->j, G->n,
                                theta, 
                                Q_store,
                                alpha,
                                0,
                                0,
                                1,
                                &ret[0], 
                                mu_store);    
        }
        break;
    }
    
  case INLA_CGENERIC_LOG_PRIOR:
    {
      ret = Calloc(1, double);

      ret[0] = 0.0;

      ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                  theta_prior_prec->x, theta);

	  break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}