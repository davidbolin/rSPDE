#include "cgeneric_defs.h"
//#include "stdio.h"
double *inla_cgeneric_rspde_intrinsic_int_model(inla_cgeneric_cmd_tp cmd, 
                                                double *theta, 
                                                inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa;
  int d, alpha;
  
  int N, M, i, k;
  
  assert(!strcasecmp(data->ints[0]->name, "n"));    
  N = data->ints[0]->ints[0];			            
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));
  int debug = data->ints[1]->ints[0];	            

  if(debug == 1){
    debug = 1;
  }

  assert(!strcasecmp(data->ints[4]->name, "d"));
  d = data->ints[4]->ints[0];
  
  assert(!strcasecmp(data->ints[5]->name, "alpha"));
  alpha = data->ints[5]->ints[0];
  
  assert(!strcasecmp(data->ints[2]->name, "graph_opt_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "graph_opt_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);

  assert(!strcasecmp(data->doubles[0]->name, "matrices_less"));
  inla_cgeneric_vec_tp *fem = data->doubles[0];
  assert(M*(alpha + 1) == fem->len);

  // prior parameters
  assert(!strcasecmp(data->doubles[1]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[1];

  assert(!strcasecmp(data->mats[0]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[0];

  assert(!strcasecmp(data->doubles[2]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[2];

  if (theta) {
        
    ltau = theta[0];
    lkappa = theta[1];
    
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

      int one = 1;
      
      if(alpha == 1){
          double sqkappatau1 = SQR(tau) * pow(kappa, 2);
          double sqkappatau2 = SQR(tau);

          dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
          dscal_(&M, &sqkappatau1, &ret[k], &one);
          daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
          
      } else if (alpha == 2){
          double sqkappatau1 = SQR(tau) * pow(kappa, 4);
          double sqkappatau2 = SQR(tau) * 2 * pow(kappa, 2);
          double sqtaukappatmp = SQR(tau);
          
          dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
          dscal_(&M, &sqkappatau1, &ret[k], &one);
          daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
          daxpy_(&M, &sqtaukappatmp, &fem->doubles[2*M], &one, &ret[k], &one);
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

     ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                theta_prior_prec->x, theta);

     //double theta_prior_mean_spde[2], theta_spde[2];
     //theta_spde[1] = lkappa;
     //theta_spde[0] = ltau;     
     //theta_prior_mean_spde[1] = 0.5 * log(8.0) - theta_prior_mean->doubles[1];
     //theta_prior_mean_spde[0] = - theta_prior_mean->doubles[0] + 0.5 *(
     //     - 2.0 * theta_prior_mean_spde[1] - (d/2.0) * log(4 * M_PI));

     //   ret[0] += logmultnormvdens(2, theta_prior_mean_spde,
     //                                 theta_prior_prec->x, theta_spde);
      
	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}