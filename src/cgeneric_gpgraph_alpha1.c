#include "cgeneric_defs.h"
#include "stdio.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_gpgraph_alpha1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  double lkappa, lsigma, kappa, sigma;

  double c1, c2, c_1, c_2, one_m_c2, l_e;

  int N, M, i, j, k;

  char *parameterization;
  
  // the size of the model
  assert(data->n_ints == 7);

  // the number of doubles (updated to include bounds and beta prior params)
  assert(data->n_doubles == 17);

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

  assert(!strcasecmp(data->ints[4]->name, "index_graph"));
  inla_cgeneric_vec_tp *idx_ij = data->ints[4];
  int M_full = idx_ij->len;

  assert(!strcasecmp(data->ints[5]->name, "count_idx"));
  inla_cgeneric_vec_tp *count_idx = data->ints[5];
  assert(M == count_idx->len);

  assert(!strcasecmp(data->ints[6]->name, "stationary_endpoints"));
  inla_cgeneric_vec_tp *stationary_endpoints = data->ints[6];

  assert(!strcasecmp(data->doubles[0]->name, "EtV2"));
  inla_cgeneric_vec_tp *EtV2 = data->doubles[0];
  
  int nE = EtV2 -> len;

  assert(!strcasecmp(data->doubles[1]->name, "EtV3"));
  inla_cgeneric_vec_tp *EtV3 = data->doubles[1];

  assert(nE == EtV3 -> len);

  assert(!strcasecmp(data->doubles[2]->name, "El"));
  inla_cgeneric_vec_tp *El = data->doubles[2];

  // prior parameters
  assert(!strcasecmp(data->doubles[3]->name, "start_theta"));
  double start_theta = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "start_lsigma"));
  double start_lsigma = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior_theta_meanlog"));
  double prior_theta_meanlog = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "prior_theta_sdlog"));
  double prior_theta_sdlog = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "prior_sigma_meanlog"));
  double prior_sigma_meanlog = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "prior_sigma_sdlog"));
  double prior_sigma_sdlog = data->doubles[8]->doubles[0];

  // Beta prior parameters for theta (kappa/range)
  assert(!strcasecmp(data->doubles[9]->name, "prior_theta_mean"));
  double prior_theta_mean = data->doubles[9]->doubles[0];

  assert(!strcasecmp(data->doubles[10]->name, "prior_theta_prec"));
  double prior_theta_prec = data->doubles[10]->doubles[0];

  // Beta prior parameters for sigma (tau/sigma)
  assert(!strcasecmp(data->doubles[11]->name, "prior_sigma_mean"));
  double prior_sigma_mean = data->doubles[11]->doubles[0];

  assert(!strcasecmp(data->doubles[12]->name, "prior_sigma_prec"));
  double prior_sigma_prec = data->doubles[12]->doubles[0];

  // Bounds for theta (kappa/range)
  assert(!strcasecmp(data->doubles[13]->name, "theta_lower_bound"));
  double theta_lower_bound = data->doubles[13]->doubles[0];

  assert(!strcasecmp(data->doubles[14]->name, "theta_upper_bound"));
  double theta_upper_bound = data->doubles[14]->doubles[0];

  // Bounds for sigma (tau/sigma)
  assert(!strcasecmp(data->doubles[15]->name, "sigma_lower_bound"));
  double sigma_lower_bound = data->doubles[15]->doubles[0];

  assert(!strcasecmp(data->doubles[16]->name, "sigma_upper_bound"));
  double sigma_upper_bound = data->doubles[16]->doubles[0];

  assert(!strcasecmp(data->chars[2]->name, "parameterization"));
  parameterization = &data->chars[2]->chars[0];

  if (theta) {
    // interpretable parameters
    double theta_val, sigma_val;

    // Transform theta (kappa/range) based on bounds
    if(theta_upper_bound > 0){
      // Upper bound exists, use logit transformation
      // theta_internal = log(x - lower) - log(upper - x)
      // => x = lower + (upper - lower) * exp(theta_internal) / (1 + exp(theta_internal))
      double exp_theta = exp(theta[1]);
      theta_val = theta_lower_bound + (theta_upper_bound - theta_lower_bound) * exp_theta / (1.0 + exp_theta);
    } else if(theta_lower_bound > 0){
      // Only lower bound, use shifted exponential
      // theta_internal = log(x - lower)
      // => x = lower + exp(theta_internal)
      theta_val = theta_lower_bound + exp(theta[1]);
    } else{
      // No bounds (or lower bound is 0), standard log transformation
      theta_val = exp(theta[1]);
    }

    // Transform sigma based on bounds
    if(sigma_upper_bound > 0){
      // Upper bound exists, use logit transformation
      double exp_sigma = exp(theta[0]);
      sigma_val = sigma_lower_bound + (sigma_upper_bound - sigma_lower_bound) * exp_sigma / (1.0 + exp_sigma);
    } else if(sigma_lower_bound > 0){
      // Only lower bound, use shifted exponential
      sigma_val = sigma_lower_bound + exp(theta[0]);
    } else{
      // No bounds (or lower bound is 0), standard log transformation
      sigma_val = exp(theta[0]);
    }

    if(!strcasecmp(parameterization, "matern")){
      // theta_val is range, sigma_val is sigma
      lkappa = log(2.0) - log(theta_val);
      kappa = exp(lkappa);
      sigma = sigma_val;
      lsigma = log(sigma);
    } else {
      // theta_val is kappa, sigma_val is tau
      kappa = theta_val;
      lkappa = log(kappa);
      // sigma_val is tau, need to convert to sigma
      lsigma = log(sigma_val);
      sigma = sigma_val;
    }
  }
  else {
    lsigma = lkappa = sigma = kappa = NAN;
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

      double *raw_entries;
      raw_entries = Calloc(M_full, double);

      // int count=0;
      //     for(i = 0; i < nE; i++){
      //       l_e = El->doubles[i];
      //       c1 = exp(-kappa*l_e);
      //       c2 = SQR(c1);
      //       one_m_c2 = 1-c2;
      //       c_1 = 0.5 + c2/one_m_c2;
      //       c_2 = -c1/one_m_c2;

      //       if(EtV2->doubles[i] != EtV3->doubles[i]){
            
      //       ret[k + idx_ij->ints[count]] = c_1;
      //       ret[k + idx_ij->ints[count + 1]] = c_1;
      //       ret[k + idx_ij->ints[count + 2]] = c_2;
      //       count += 3;
      //     }else{
      //       ret[k + idx_ij->ints[count]] = tanh(0.5 * kappa * l_e);
      //       count++;
      //     }
      //   }

      //   if(stationary_endpoints->ints[0]>=0){
      //     int stat_endpt_len = stationary_endpoints->len;
      //     for(i = 0; i < stat_endpt_len; i++){
      //       ret[k+idx_ij->ints[count+i]] = 0.5;
      //     }
      //     count += stat_endpt_len;
      //   }

      //   assert(count == M);

      //   double fact = 2*kappa / (pow(sigma,2));

      //   int one=1;
      //   dscal_(&M, &fact, &ret[k], &one);
    
       int count=0;
          for(i = 0; i < nE; i++){
            l_e = El->doubles[i];
            c1 = exp(-kappa*l_e);
            c2 = SQR(c1);
            one_m_c2 = 1-c2;
            c_1 = 0.5 + c2/one_m_c2;
            c_2 = -c1/one_m_c2;

            if(EtV2->doubles[i] != EtV3->doubles[i]){
            
            raw_entries[idx_ij->ints[count]] = c_1;
            raw_entries[idx_ij->ints[count + 1]] = c_1;
            raw_entries[idx_ij->ints[count + 2]] = c_2;
            count += 3;
          }else{
            raw_entries[idx_ij->ints[count]] = tanh(0.5 * kappa * l_e);
            count++;
          }
        }

        if(stationary_endpoints->ints[0]>=0){
          int stat_endpt_len = stationary_endpoints->len;
          for(i = 0; i < stat_endpt_len; i++){
            raw_entries[idx_ij->ints[count+i]] = 0.5;
          }
          count += stat_endpt_len;
        }

        assert(count == M_full);

        double fact = 2*kappa / (pow(sigma,2));

        int one=1;
        dscal_(&M_full, &fact, &raw_entries[0], &one);

        count = 0;
        for(i = 0; i < M; i++){
          for(j = 0; j < count_idx->ints[i]; j++){
            // ret[k + i] += raw_entries[count];
            ret[k + i] += raw_entries[count];
            count++;
          }
        }
        assert(M_full == count);

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

      // Prior for theta (kappa/range)
      if(theta_upper_bound > 0){
        // Beta prior for bounded parameter
        double theta_val;
        double exp_theta = exp(theta[1]);
        theta_val = theta_lower_bound + (theta_upper_bound - theta_lower_bound) * exp_theta / (1.0 + exp_theta);

        double s_1 = (prior_theta_mean - theta_lower_bound) / (theta_upper_bound - theta_lower_bound) * prior_theta_prec;
        double s_2 = (1 - (prior_theta_mean - theta_lower_bound) / (theta_upper_bound - theta_lower_bound)) * prior_theta_prec;
        double x_scaled = (theta_val - theta_lower_bound) / (theta_upper_bound - theta_lower_bound);

        ret[0] += logdbeta(x_scaled, s_1, s_2) - log(theta_upper_bound - theta_lower_bound);
        // Add Jacobian adjustment for logit transformation
        ret[0] += theta[1] - 2.0 * log(1.0 + exp_theta);
      } else if(theta_lower_bound > 0){
        // Shifted exponential, use normal prior on log scale
        ret[0] += -0.5 * SQR(theta[1] - prior_theta_meanlog)/(SQR(prior_theta_sdlog)) -
          log(prior_theta_sdlog) - 0.5 * log(2.0 * M_PI);
      } else{
        // Standard log transformation, use normal prior
        ret[0] += -0.5 * SQR(theta[1] - prior_theta_meanlog)/(SQR(prior_theta_sdlog)) -
          log(prior_theta_sdlog) - 0.5 * log(2.0 * M_PI);
      }

      // Prior for sigma (tau/sigma)
      if(sigma_upper_bound > 0){
        // Beta prior for bounded parameter
        double sigma_val;
        double exp_sigma = exp(theta[0]);
        sigma_val = sigma_lower_bound + (sigma_upper_bound - sigma_lower_bound) * exp_sigma / (1.0 + exp_sigma);

        double s_1 = (prior_sigma_mean - sigma_lower_bound) / (sigma_upper_bound - sigma_lower_bound) * prior_sigma_prec;
        double s_2 = (1 - (prior_sigma_mean - sigma_lower_bound) / (sigma_upper_bound - sigma_lower_bound)) * prior_sigma_prec;
        double x_scaled = (sigma_val - sigma_lower_bound) / (sigma_upper_bound - sigma_lower_bound);

        ret[0] += logdbeta(x_scaled, s_1, s_2) - log(sigma_upper_bound - sigma_lower_bound);
        // Add Jacobian adjustment for logit transformation
        ret[0] += theta[0] - 2.0 * log(1.0 + exp_sigma);
      } else if(sigma_lower_bound > 0){
        // Shifted exponential, use normal prior on log scale
        ret[0] += -0.5 * SQR(theta[0] - prior_sigma_meanlog)/(SQR(prior_sigma_sdlog)) -
          log(prior_sigma_sdlog) - 0.5 * log(2.0 * M_PI);
      } else{
        // Standard log transformation, use normal prior
        ret[0] += -0.5 * SQR(lsigma - prior_sigma_meanlog)/(SQR(prior_sigma_sdlog)) -
          log(prior_sigma_sdlog) - 0.5 * log(2.0 * M_PI);
      }
	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}