#include "cgeneric_defs.h"
// #include "stdio.h"
// #include "gsl/gsl_vector_double.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_stat_frac_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa;
  double alpha, nu;
  int m_alpha;
  int N, M, i, k, j, rspde_order, d;
  int full_size, less_size;
  int one = 1;
  char *parameterization, *theta_param;

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

  assert(!strcasecmp(data->ints[4]->name, "rspde.order"));
  rspde_order = data->ints[4]->ints[0];

  assert(!strcasecmp(data->ints[5]->name, "d"));
  d = data->ints[5]->ints[0];

  assert(!strcasecmp(data->chars[2]->name, "parameterization"));
  parameterization = &data->chars[2]->chars[0];

  assert(!strcasecmp(data->chars[3]->name, "prior.theta.param"));
  theta_param = &data->chars[3]->chars[0];

  assert(!strcasecmp(data->doubles[0]->name, "nu"));
  nu = data->doubles[0]->doubles[0];

  alpha = nu + d / 2.0;
  m_alpha = floor(alpha);

  assert(!strcasecmp(data->doubles[1]->name, "matrices_less"));
  inla_cgeneric_vec_tp *fem_less = data->doubles[1];

  assert(!strcasecmp(data->doubles[2]->name, "matrices_full"));
  inla_cgeneric_vec_tp *fem_full = data->doubles[2];
  full_size = (fem_full->len)/(m_alpha+2);
  less_size = (fem_less->len)/(m_alpha+1);
  assert(M == rspde_order * full_size + less_size);


  assert(!strcasecmp(data->doubles[3]->name, "r_ratapprox"));
  double *r = data->doubles[3]->doubles;

  assert(!strcasecmp(data->doubles[4]->name, "p_ratapprox"));
  double *p = data->doubles[4]->doubles;
  
  assert(!strcasecmp(data->doubles[5]->name, "k_ratapprox"));
  double k_rat = data->doubles[5]->doubles[0];

  // prior parameters
  
  assert(!strcasecmp(data->doubles[6]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[6];

  assert(!strcasecmp(data->mats[0]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[0];

  assert(!strcasecmp(data->doubles[7]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[7];

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

      // FORTRAN IMPLEMENTATION

      double multQ = pow(kappa, 2*alpha) * SQR(tau);

      switch(m_alpha){
        case 0:
        {
          double fact_mult;
          for(j = 0; j < rspde_order; j++){
            fact_mult = multQ * (1-p[j])/r[j];
            dcopy_(&full_size, &fem_full->doubles[0], &one, &ret[k + j*full_size], &one);
            dscal_(&full_size, &fact_mult, &ret[k+ j*full_size], &one);
            fact_mult = multQ / (r[j] * SQR(kappa));
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[full_size], &one, &ret[k+j*full_size], &one);
          }
          // dcopy_(&less_size, &fem_less->doubles[0], &one, &ret[k+rspde_order * full_size], &one);
          // fact_mult = multQ/k_rat;
          // dscal_(&less_size, &fact_mult, &ret[k+rspde_order * full_size], &one);
            for(i = 0; i < less_size; i++){
              if(fem_less->doubles[i] != 0){
                ret[k+rspde_order*full_size + i] = multQ * (
                    1/(k_rat * fem_less->doubles[i])
                );
              }
            }          
          break;
        }
        case 1:
        {
          double *Malpha2, fact_mult;
          Malpha2 = Calloc(full_size, double);
          for(j = 0; j < rspde_order; j++){
            dcopy_(&full_size, &fem_full->doubles[0], &one, &ret[k+j*full_size], &one);
            fact_mult = 1/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[full_size], &one, &ret[k+j*full_size], &one);
            fact_mult = multQ * (1-p[j])/r[j];
            dscal_(&full_size, &fact_mult, &ret[k+j*full_size], &one);
            dcopy_(&full_size, &fem_full->doubles[full_size], &one, Malpha2, &one);
            fact_mult = 1/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[2*full_size], &one, Malpha2, &one);
            fact_mult = multQ/(SQR(kappa) * r[j]);
            daxpy_(&full_size, &fact_mult, Malpha2, &one, &ret[k + j*full_size], &one);
          }

          free(Malpha2);

          dcopy_(&less_size, &fem_less->doubles[0], &one, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ/k_rat;
          dscal_(&less_size, &fact_mult, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ/(k_rat * SQR(kappa));
          daxpy_(&less_size, &fact_mult, &fem_less->doubles[less_size], &one, &ret[k+rspde_order*full_size], &one);
          break;
        }
        default:
        {
          double *Malpha2, fact_mult;
          Malpha2 = Calloc(full_size, double);
          for(j = 0; j < rspde_order; j++){
            dcopy_(&full_size, &fem_full->doubles[0],&one, &ret[k+j*full_size], &one);
            fact_mult = m_alpha/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[full_size], &one, &ret[k+j*full_size], &one);
            for(i = 2; i<= m_alpha; i++){
              fact_mult = nChoosek(m_alpha, i)/(pow(kappa, 2*i));
              daxpy_(&full_size, &fact_mult, &fem_full->doubles[i*full_size], &one, &ret[k+j*full_size], &one);
            }
            fact_mult = multQ * (1-p[j])/r[j];
            dscal_(&full_size, &fact_mult, &ret[k+j*full_size], &one);
            dcopy_(&full_size, &fem_full->doubles[full_size], &one, Malpha2, &one);
            fact_mult = m_alpha/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[2*full_size], &one, Malpha2, &one);
            for(i = 2; i<= m_alpha; i++){
              fact_mult = nChoosek(m_alpha, i)/(pow(kappa, 2*i));
              daxpy_(&full_size, &fact_mult, &fem_full->doubles[(i+1)*full_size], &one, Malpha2, &one);
            }
            fact_mult = multQ/(SQR(kappa) * r[j]);
            daxpy_(&full_size, &fact_mult, Malpha2, &one, &ret[k + j*full_size], &one);
          }

          free(Malpha2);

          dcopy_(&less_size, &fem_less->doubles[0], &one, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ/k_rat;
          dscal_(&less_size, &fact_mult, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ * m_alpha/(k_rat * SQR(kappa));
          daxpy_(&less_size, &fact_mult, &fem_less->doubles[less_size], &one, &ret[k+rspde_order*full_size], &one);
          for(j = 2; j<= m_alpha; j++){
            fact_mult = multQ * nChoosek(m_alpha, j)/(k_rat * pow(SQR(kappa),j));
            daxpy_(&less_size, &fact_mult, &fem_less->doubles[j*less_size], &one, &ret[k+rspde_order*full_size], &one);            
          }
          break;
        }
      }

      // GSL IMPLEMENTATION

      // gsl_vector * FEM1 = gsl_vector_calloc(full_size); // C, then G, G_2, etc.
      // gsl_vector * FEM2 = gsl_vector_calloc(full_size); // G, then G_2, G_3, etc., then part to be returned



      // switch(new_m_alpha){
      //   case 0:
      //   {
      //     dcopy_(&full_size, &fem_full->doubles[0], &one, &FEM1->data[0], &one);
      //     for(j = 0; j < rspde_order; j++){
      //       dcopy_(&full_size, &fem_full->doubles[full_size], &one, &FEM2->data[0], &one);
      //       gsl_vector_axpby(multQ * (1-p[j])/r[j], FEM1, multQ / (r[j] * SQR(kappa)), FEM2);
      //       dcopy_(&full_size, &FEM2->data[0], &one, &ret[k + j*full_size], &one);
      //     }
      //     gsl_vector_free(FEM1);
      //     gsl_vector_free(FEM2);
      //     gsl_vector * FEM1 = gsl_vector_calloc(less_size); // C, then G, G_2, etc.
      //     dcopy_(&less_size, &fem_less->doubles[0], &one, &FEM1->data[0], &one);
      //     gsl_vector_scale(FEM1, multQ/k_rat);
      //     dcopy_(&less_size, &FEM1->data[0], &one, &ret[k + rspde_order * full_size], &one);
      //     gsl_vector_free(FEM1);
      //     break;
      //   }
      //   case 1:
      //   {
      //     dcopy_(&full_size, &fem_full->doubles[0], &one, &FEM1->data[0], &one);
      //     dcopy_(&full_size, &fem_full->doubles[full_size], &one, &FEM2->data[0], &one);
      //     gsl_vector_axpby(1/SQR(kappa), FEM2, 1, FEM1); // FEM1 = M_alpha
      //     gsl_vector * FEM3 = gsl_vector_calloc(full_size); 
      //     dcopy_(&full_size, &fem_full->doubles[2 * full_size], &one, &FEM3->data[0], &one);
      //     gsl_vector_axpby(1/SQR(kappa), FEM3, 1, FEM2); // FEM2 = M_alpha2
      //     gsl_vector_memcpy(FEM3, FEM2);
      //     for(j = 0; j < rspde_order; j++){
      //           gsl_vector_axpby(multQ * (1-p[j])/r[j], FEM1, multQ/(SQR(kappa) * r[j]), FEM2); //FEM2 -> part of Q
      //           dcopy_(&full_size, &FEM2->data[0], &one, &ret[k + j*full_size], &one);
      //           gsl_vector_memcpy(FEM2, FEM3);
      //     }
      //     // Add k part
      //     gsl_vector_free(FEM1);
      //     gsl_vector_free(FEM2);
      //     gsl_vector_free(FEM3);
      //     gsl_vector * FEM1 = gsl_vector_calloc(less_size);
      //     gsl_vector * FEM2 = gsl_vector_calloc(less_size);
      //     dcopy_(&less_size, &fem_less->doubles[0], &one, &FEM1->data[0], &one); // copy C back to FEM1
      //     dcopy_(&less_size, &fem_less->doubles[less_size], &one, &FEM2->data[0], &one);
      //     gsl_vector_axpby(multQ/(k_rat * SQR(kappa)), FEM2, multQ/k_rat, FEM1);
      //     dcopy_(&less_size, &FEM1->data[0], &one, &ret[k + rspde_order * full_size], &one);
      //     gsl_vector_free(FEM1);
      //     gsl_vector_free(FEM2);
      //     break;
      //   }
      //   default:
      //   {
      //     dcopy_(&full_size, &fem_full->doubles[0], &one, &FEM1->data[0], &one);
      //     gsl_vector * FEM3 = gsl_vector_calloc(full_size); 
      //     dcopy_(&full_size, &fem_full->doubles[full_size], &one, &FEM2->data[0], &one);
      //     gsl_vector_axpby(new_m_alpha/SQR(kappa), FEM2, 1, FEM1); // FEM1 = M_alpha
      //     for(j = 2; j<= new_m_alpha; j++){
      //       dcopy_(&full_size, &fem_full->doubles[j * full_size], &one, &FEM3->data[0], &one);
      //       gsl_vector_axpby(nChoosek(new_m_alpha, j)/(pow(kappa, 2*j)), FEM3, 1, FEM1);
      //     }
      //     dcopy_(&full_size, &fem_full->doubles[2 * full_size], &one, &FEM3->data[0], &one);
      //     gsl_vector_axpby(new_m_alpha/SQR(kappa), FEM3, 1, FEM2); // FEM2 = M_alpha2
      //     for(j = 2; j<= new_m_alpha; j++){
      //       dcopy_(&full_size, &fem_full->doubles[(j+1) * full_size], &one, &FEM3->data[0], &one);
      //       gsl_vector_axpby(nChoosek(new_m_alpha, j)/(pow(kappa,2*j)), FEM3, 1, FEM2);
      //     }


      //     gsl_vector_memcpy(FEM3, FEM2);
      //     for(j = 0; j < rspde_order; j++){
      //           gsl_vector_axpby(multQ * (1-p[j])/r[j], FEM1, multQ/(SQR(kappa) * r[j]), FEM2); //FEM2 -> part of Q
      //           dcopy_(&full_size, &FEM2->data[0], &one, &ret[k + j*full_size], &one);
      //           gsl_vector_memcpy(FEM2, FEM3);
      //     }
      //     // Add k part
      //     gsl_vector_free(FEM1);
      //     gsl_vector_free(FEM2);
      //     gsl_vector_free(FEM3);
      //     gsl_vector * FEM1 = gsl_vector_calloc(less_size);
      //     gsl_vector * FEM2 = gsl_vector_calloc(less_size);
      //     dcopy_(&less_size, &fem_less->doubles[0], &one, &FEM1->data[0], &one); // copy C back to FEM1
      //     dcopy_(&less_size, &fem_less->doubles[less_size], &one, &FEM2->data[0], &one);
      //     gsl_vector_axpby(multQ * new_m_alpha/(k_rat * SQR(kappa)), FEM2, multQ/k_rat, FEM1);
      //     for(j = 2; j<= new_m_alpha; j++){
      //       dcopy_(&less_size, &fem_less->doubles[j * less_size], &one, &FEM2->data[0], &one);
      //       gsl_vector_axpby(multQ * nChoosek(new_m_alpha, j)/(k_rat * pow(SQR(kappa),j)), FEM2, 1, FEM1);
      //     }
      //     dcopy_(&less_size, &FEM1->data[0], &one, &ret[k + rspde_order * full_size], &one);
      //     gsl_vector_free(FEM1);
      //     gsl_vector_free(FEM2);
      //     break;
      //   }
      // }

      // DIRECT C IMPLEMENTATION

      // double *Malpha, *Malpha2;


      // // double multQ = pow(kappa, 2*new_alpha) * SQR(tau);

      // if(new_m_alpha == 0){
      //   for(j = 0; j < rspde_order; j++){
      //       for(i = 0; i < full_size; i++){
      //           ret[k + j*full_size + i] = multQ * (
      //               (fem_full->doubles[i] + (fem_full->doubles[full_size+i])/(SQR(kappa))) -
      //                (p[j] * fem_full->doubles[i])
      //           ) / r[j];
      //       }
      //   }

      //   // Kpart
      //   for(i = 0; i < less_size; i++){
      //           ret[k+rspde_order*full_size + i] = multQ * (
      //               fem_less->doubles[i]/k_rat
      //           );
      //       }

      // } else{

      // Malpha = Calloc(full_size, double);
      // Malpha2 = Calloc(full_size, double);

      // if(new_m_alpha == 1){
      //   // for(i = 0; i < full_size; i++){
      //   //     Malpha[i] = fem_full->doubles[i] + (fem_full->doubles[full_size+i])/(SQR(kappa));
      //   //     Malpha2[i] = fem_full->doubles[full_size+i] + (fem_full->doubles[2*full_size+i])/(SQR(kappa));
      //   // }
      //   for(i = 0; i < full_size; i++){
      //       Malpha[i] = fem_full->doubles[i] + (fem_full->doubles[full_size+i])/(SQR(kappa));
      //   }
      //   for(i = 0; i < full_size; i++){
      //       Malpha2[i] = fem_full->doubles[full_size+i] + (fem_full->doubles[2*full_size+i])/(SQR(kappa));
      //   }
      // } else if(new_m_alpha > 1){
      //   // for(i = 0; i < full_size; i++){
      //   //     Malpha[i] = fem_full->doubles[i] + new_m_alpha * (fem_full->doubles[full_size+i])/(SQR(kappa));
      //   //     for(j = 2; j <= new_m_alpha; j++){
      //   //         Malpha[i] += nChoosek(new_m_alpha, j) * (fem_full->doubles[j*full_size+i])/(pow(kappa, 2*j));
      //   //     }

      //   //     Malpha2[i] = fem_full->doubles[full_size+i] + new_m_alpha * (fem_full->doubles[2*full_size+i])/(SQR(kappa));
      //   //     for(j = 2; j <= new_m_alpha ; j++){
      //   //         Malpha2[i] += nChoosek(new_m_alpha, j) * (fem_full->doubles[(j+1)*full_size + i])/(pow(kappa,2*j));
      //   //     }
      //   // }

      //   for(i = 0; i < full_size; i++){
      //       Malpha[i] = fem_full->doubles[i] + new_m_alpha * (fem_full->doubles[full_size+i])/(SQR(kappa));
      //       for(j = 2; j <= new_m_alpha; j++){
      //           Malpha[i] += nChoosek(new_m_alpha, j) * (fem_full->doubles[j*full_size+i])/(pow(kappa, 2*j));
      //       }
      //   }
      //   for(i = 0; i < full_size; i++){
      //       Malpha2[i] = fem_full->doubles[full_size+i] + new_m_alpha * (fem_full->doubles[2*full_size+i])/(SQR(kappa));
      //       for(j = 2; j <= new_m_alpha ; j++){
      //           Malpha2[i] += nChoosek(new_m_alpha, j) * (fem_full->doubles[(j+1)*full_size + i])/(pow(kappa,2*j));
      //       }
      //   }
      // }

      // for(j = 0; j < rspde_order; j++){
      //   for(i = 0; i < full_size; i++){
      //       ret[k + j * full_size + i] = multQ * (
      //           (1-p[j]) * Malpha[i] + (Malpha2[i])/(SQR(kappa))
      //       )/(r[j]);
      //   }
      // }

      // if(new_m_alpha == 1){
        
      //   for(i = 0; i < less_size; i++){
      //   ret[k+rspde_order*full_size + i] = multQ/(k_rat) * (
      //               fem_less->doubles[i] + (fem_less->doubles[less_size+i])/(SQR(kappa))
      //   );
      //   }
      // } else{

      //   for(i = 0; i < less_size; i++){
      //       ret[k+rspde_order*full_size + i] = multQ/(k_rat) * (
      //           fem_less->doubles[i] + (new_m_alpha/SQR(kappa)) * fem_less->doubles[less_size+i]
      //       );
      //       for(j = 2; j <= new_m_alpha ; j++){
      //       ret[k+rspde_order*full_size + i] += multQ/(k_rat) * (
      //                   nChoosek(new_m_alpha,j)*(fem_less->doubles[i + j*less_size])/(pow(kappa,2*j))
      //       );
      //       }
      //       }
      //   }

      // }

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