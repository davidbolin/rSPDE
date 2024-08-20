#include "cgeneric_defs.h"
#include "stdio.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_gpgraph_alpha2_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  double lkappa, lsigma, kappa, sigma, tau;

  int N, M, k, i;

  char *parameterization;
  
  // the size of the model
  assert(data->n_ints == 8);

  // the number of doubles
  assert(data->n_doubles == 7);

  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug == 1){
    debug = 1;
  }

  assert(!strcasecmp(data->ints[2]->name, "prec_graph_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "prec_graph_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);

  assert(!strcasecmp(data->ints[4]->name, "upper_edges"));
  inla_cgeneric_vec_tp *upper_edges = data->ints[4];

  assert(!strcasecmp(data->ints[5]->name, "lower_edges"));
  inla_cgeneric_vec_tp *lower_edges = data->ints[5];

  assert(!strcasecmp(data->ints[6]->name, "upper_edges_len"));
  int upper_edges_len = data->ints[6]->ints[0];

  assert(!strcasecmp(data->ints[7]->name, "lower_edges_len"));
  int lower_edges_len = data->ints[7]->ints[0];

  assert(!strcasecmp(data->smats[0]->name, "Tc"));
  inla_cgeneric_smat_tp *Tc = data->smats[0];

  assert(!strcasecmp(data->doubles[0]->name, "El"));
  inla_cgeneric_vec_tp *El = data->doubles[0];  
  
  int nE = El -> len;

  // prior parameters
  assert(!strcasecmp(data->doubles[1]->name, "start_theta"));
  double start_theta = data->doubles[1]->doubles[0];

  assert(!strcasecmp(data->doubles[2]->name, "start_lsigma"));
  double start_lsigma = data->doubles[2]->doubles[0];

  assert(!strcasecmp(data->doubles[3]->name, "prior_theta_meanlog"));
  double prior_theta_meanlog = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "prior_theta_sdlog"));
  double prior_theta_sdlog = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior_sigma_meanlog"));
  double prior_sigma_meanlog = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "prior_sigma_sdlog"));
  double prior_sigma_sdlog = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->chars[2]->name, "parameterization"));
  parameterization = &data->chars[2]->chars[0];

  if (theta) {
    // interpretable parameters 

    if(!strcasecmp(parameterization, "matern")){
      lkappa =  0.5 * log(12.0) - theta[1];
    } else {
      lkappa = theta[1];
    }
    lsigma = theta[0];
    kappa = exp(lkappa);
    sigma = exp(lsigma);
    tau = 1.0/sigma;
  }
  else {   
    lsigma = lkappa = sigma = kappa = tau = NAN;
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
      ret[0] = N;       /* dimension */
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
      
      compute_Q_alpha2(Tc->i, Tc->j, Tc->x, kappa, tau, nE, 0.5,
                            Tc->nrow, Tc->ncol, Tc->n, El->doubles, &ret[k], lower_edges->ints,
                                        upper_edges->ints, lower_edges_len, upper_edges_len);

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
      ret[1] = start_lsigma;
      ret[2] = start_theta;
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

      ret[0] += -0.5 * SQR(theta[1] - prior_theta_meanlog)/(SQR(prior_theta_sdlog)) - 
      log(prior_theta_sdlog) - 0.5 * log(2.0 * M_PI); 

      ret[0] += -0.5 * SQR(lsigma - prior_sigma_meanlog)/(SQR(prior_sigma_sdlog)) - 
      log(prior_sigma_sdlog) - 0.5 * log(2.0 * M_PI);
	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}