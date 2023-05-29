#include "cgeneric_defs.h"
#include "stdio.h"
// #include "gsl/gsl_vector_double.h"

double cut_decimals(double nu){
    double temp = nu - floor(nu);
    if(temp < pow(10,-3)){
        temp = pow(10,-3);
    }
    if(temp > 0.999){
        temp = 0.999;
    }
    return temp;
}

double pnorm(double x, double mu, double sd) 
{
    return (1 + erf((x-mu) / (sd * sqrt(2.0))))/(2.0);
}

double logdbeta(double x, double s_1, double s_2){
    double tmp = lgamma(s_1 + s_2) - lgamma(s_1) - lgamma(s_2);
    tmp += (s_1-1)*log(x) + (s_2-1)*log(1-x);
    return tmp;
}

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_stat_general_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa, lnu, nu;
  double alpha, nu_upper_bound;
  int m_alpha;
  double prior_nu_mean, prior_nu_loglocation, prior_nu_prec;
  double prior_nu_logscale;
  double start_nu;
  int N, M, i, k, j, rspde_order;
  double d;
  char *prior_nu_dist, *parameterization, *theta_param;
  int full_size, less_size;
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

  assert(!strcasecmp(data->ints[4]->name, "rspde.order"));
  rspde_order = data->ints[4]->ints[0];

  assert(!strcasecmp(data->chars[2]->name, "prior.nu.dist"));
  prior_nu_dist = &data->chars[2]->chars[0];

  assert(!strcasecmp(data->chars[4]->name, "prior.theta.param"));
  theta_param = &data->chars[4]->chars[0];

  assert(!strcasecmp(data->chars[3]->name, "parameterization"));
  parameterization = &data->chars[3]->chars[0];

  assert(!strcasecmp(data->doubles[0]->name, "d"));
  d = data->doubles[0]->doubles[0];

  assert(!strcasecmp(data->doubles[1]->name, "nu.upper.bound"));
  nu_upper_bound = data->doubles[1]->doubles[0];

  alpha = nu_upper_bound + d / 2.0;
  m_alpha = floor(alpha);

  assert(!strcasecmp(data->doubles[2]->name, "matrices_less"));
  inla_cgeneric_vec_tp *fem_less = data->doubles[2];

  assert(!strcasecmp(data->doubles[3]->name, "matrices_full"));
  inla_cgeneric_vec_tp *fem_full = data->doubles[3];
  full_size = (fem_full->len)/(m_alpha+2);
  less_size = (fem_less->len)/(m_alpha+1);
  assert(M == rspde_order * full_size + less_size);


  assert(!strcasecmp(data->mats[0]->name, "rational_table"));
  inla_cgeneric_mat_tp *rational_table = data->mats[0];
  assert(rational_table->nrow == 999);  

  // prior parameters
  assert(!strcasecmp(data->doubles[4]->name, "start.theta"));
  inla_cgeneric_vec_tp *start_theta = data->doubles[4];
  
  assert(!strcasecmp(data->doubles[5]->name, "theta.prior.mean"));
  inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[5];

  assert(!strcasecmp(data->mats[1]->name, "theta.prior.prec"));
  inla_cgeneric_mat_tp *theta_prior_prec = data->mats[1];

  assert(!strcasecmp(data->doubles[6]->name, "prior.nu.loglocation"));
  prior_nu_loglocation = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "prior.nu.mean"));
  prior_nu_mean = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "prior.nu.prec"));
  prior_nu_prec = data->doubles[8]->doubles[0];

  assert(!strcasecmp(data->doubles[9]->name, "prior.nu.logscale"));
  prior_nu_logscale = data->doubles[9]->doubles[0];

  assert(!strcasecmp(data->doubles[10]->name, "start.nu"));
  start_nu = data->doubles[10]->doubles[0];

  if (theta) {
    // interpretable parameters 
    lnu = theta[2];
    nu = (exp(lnu)/(1.0 + exp(lnu))) * nu_upper_bound;
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
    ltau = lkappa = lnu = tau = kappa = nu = NAN;
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

      double new_alpha = nu + d / 2.0;

      int new_m_alpha = (int) floor(new_alpha);

      double multQ = pow(kappa, 2*new_alpha) * SQR(tau);

      int row_nu = (int)round(1000*cut_decimals(new_alpha))-1;

      double *rat_coef = Calloc(n_terms-1, double);
      
      rat_coef = &rational_table->x[row_nu*n_terms+1];

      double *r, *p, k_rat;

      r = Calloc(rspde_order, double);
      p = Calloc(rspde_order, double);
      
      r = &rat_coef[0];
      p = &rat_coef[rspde_order];
      k_rat = rat_coef[2*rspde_order];

      // FORTRAN IMPLEMENTATION



      switch(new_m_alpha){
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
            fact_mult = new_m_alpha/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[full_size], &one, &ret[k+j*full_size], &one);
            for(i = 2; i<= new_m_alpha; i++){
              fact_mult = nChoosek(new_m_alpha, i)/(pow(kappa, 2*i));
              daxpy_(&full_size, &fact_mult, &fem_full->doubles[i*full_size], &one, &ret[k+j*full_size], &one);
            }
            fact_mult = multQ * (1-p[j])/r[j];
            dscal_(&full_size, &fact_mult, &ret[k+j*full_size], &one);
            dcopy_(&full_size, &fem_full->doubles[full_size], &one, Malpha2, &one);
            fact_mult = new_m_alpha/SQR(kappa);
            daxpy_(&full_size, &fact_mult, &fem_full->doubles[2*full_size], &one, Malpha2, &one);
            for(i = 2; i<= new_m_alpha; i++){
              fact_mult = nChoosek(new_m_alpha, i)/(pow(kappa, 2*i));
              daxpy_(&full_size, &fact_mult, &fem_full->doubles[(i+1)*full_size], &one, Malpha2, &one);
            }
            fact_mult = multQ/(SQR(kappa) * r[j]);
            daxpy_(&full_size, &fact_mult, Malpha2, &one, &ret[k + j*full_size], &one);
          }

          free(Malpha2);

          dcopy_(&less_size, &fem_less->doubles[0], &one, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ/k_rat;
          dscal_(&less_size, &fact_mult, &ret[k+rspde_order*full_size], &one);
          fact_mult = multQ * new_m_alpha/(k_rat * SQR(kappa));
          daxpy_(&less_size, &fact_mult, &fem_less->doubles[less_size], &one, &ret[k+rspde_order*full_size], &one);
          for(j = 2; j<= new_m_alpha; j++){
            fact_mult = multQ * nChoosek(new_m_alpha, j)/(k_rat * pow(SQR(kappa),j));
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
      ret = Calloc(4, double);
      ret[0] = 3;
      ret[1] = start_theta->doubles[0];
      ret[2] = start_theta->doubles[1];
      ret[3] = log(start_nu/(nu_upper_bound - start_nu));
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

      if(!strcasecmp(theta_param, "theta") || !strcasecmp(parameterization, "spde")){
          ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                      theta_prior_prec->x, theta);
      }
      else {
        double theta_prior_mean_spde[2], theta_spde[2], prior_nu_tmp;
        if(!strcasecmp(prior_nu_dist, "lognormal")){
          prior_nu_tmp = exp(prior_nu_loglocation);
        }
        else{
          prior_nu_tmp = prior_nu_mean;
        }
        theta_spde[1] = lkappa;
        theta_spde[0] = ltau;
        theta_prior_mean_spde[1] = 0.5 * log(8.0 * prior_nu_tmp) - theta_prior_mean->doubles[1];
        theta_prior_mean_spde[0] = - theta_prior_mean->doubles[0] + 0.5 *(
          lgamma(prior_nu_tmp) - 2.0 * prior_nu_tmp * theta_prior_mean_spde[1] - 
          (d/2.0) * log(4 * M_PI) - lgamma(prior_nu_tmp + d/2.0)
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