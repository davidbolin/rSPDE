#include "cgeneric_defs.h"
#include "stdio.h"
// #include "gsl/gsl_vector_double.h"

unsigned nChoosek( int n, int k ){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_stat_int_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;
  double ltau, lkappa, tau, kappa, prior_kappa_meanlog;
  double prior_kappa_sdlog, prior_tau_meanlog, prior_tau_sdlog;
  double start_ltau, start_lkappa;

  int N, M, i, k, j;
  
  // the size of the model
  assert(data->n_ints == 5);

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

  assert(!strcasecmp(data->ints[2]->name, "m_alpha"));
  int m_alpha = data->ints[2]->ints[0];

  assert(!strcasecmp(data->ints[3]->name, "graph_opt_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[3];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[4]->name, "graph_opt_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[4];
  assert(M == graph_j->len);

  // assert(!strcasecmp(data->ints[5]->name, "positions_C"));
  // inla_cgeneric_vec_tp *positions_C = data->ints[5];

  // assert(!strcasecmp(data->ints[6]->name, "positions_G"));
  // inla_cgeneric_vec_tp *positions_G = data->ints[6];

  assert(!strcasecmp(data->doubles[0]->name, "matrices_less"));
  inla_cgeneric_vec_tp *fem = data->doubles[0];
  assert(M*(m_alpha+1) == fem->len);

  // prior parameters
  assert(!strcasecmp(data->doubles[1]->name, "prior.kappa.meanlog"));
  prior_kappa_meanlog = data->doubles[1]->doubles[0];

  assert(!strcasecmp(data->doubles[2]->name, "prior.kappa.sdlog"));
  prior_kappa_sdlog = data->doubles[2]->doubles[0];

  assert(!strcasecmp(data->doubles[3]->name, "prior.tau.meanlog"));
  prior_tau_meanlog = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "prior.tau.sdlog"));
  prior_tau_sdlog = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "start.lkappa"));
  start_lkappa = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "start.ltau"));
  start_ltau = data->doubles[6]->doubles[0];

  if (theta) {
    // interpretable parameters 
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

      // int one = 1;
      // gsl_vector * GC = gsl_vector_calloc (M); // First G then C
      // gsl_vector * retG2 = gsl_vector_calloc (M); // first G2 then return

      // dcopy_(&M, &fem->doubles[2*M], &one, &retG2->data[0], &one);
      // dcopy_(&M, &fem->doubles[M], &one, &GC->data[0], &one);
      // gsl_vector_axpby(2*SQR(tau)*SQR(kappa), GC, SQR(tau), retG2);
      // dcopy_(&M, &fem->doubles[0], &one, &GC->data[0], &one);
      // gsl_vector_axpby(SQR(tau)*SQR(kappa*kappa), GC, 1, retG2);
      // dcopy_(&M, &retG2->data[0], &one, &ret[k], &one);

      
      // gsl_vector * retV = gsl_vector_calloc(M);
      // gsl_vector_memcpy(retV, G2);

      // gsl_vector_axpby(2*SQR(tau)*SQR(kappa), G, SQR(tau), retV);
      // gsl_vector_axpby(SQR(tau)*SQR(kappa*kappa), C, 1, retV);
      // dcopy_(&M, &retV->data[0], &one, &ret[k], &one);


      //   dscal_(&M, &sqtau, &ret[k], &one); 
      //   for(i = 0; i < positions_C->len; i++){
      //     ret[k + positions_G->ints[i]-1] += 2*SQR(tau) * SQR(kappa) * fem->doubles[M+positions_G->ints[i]-1];
      //     ret[k + positions_C->ints[i]-1] += SQR(tau) * SQR(kappa * kappa) * fem->doubles[positions_C->ints[i]-1];
      //   }
      //   for(i = positions_C->len; i<positions_G->len; i++){
      //     ret[k + positions_G->ints[i]-1] += 2*SQR(tau) * SQR(kappa) * fem->doubles[M+positions_G->ints[i]-1];
      //   }

      // Direct version:

      // if(m_alpha == 1){
      //   for (i = 0; i < M; i++) {
      //     ret[k + i] = SQR(tau) * (SQR(kappa) * fem->doubles[i] + fem->doubles[M+i]);
      //   }
      // }
      // else if(m_alpha > 1){
      //   for (i = 0; i < M; i++) {
      //     ret[k + i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[i] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[M+i]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(j+2)*M+i]);
      //       }
      //     }

      //   }
      // }

      // Currently the faster version:

      // if(m_alpha == 1){
      //   for (i = 0; i < M; i++) {
      //     ret[k + i] = SQR(tau) * (SQR(kappa) * fem->doubles[2*i] + fem->doubles[2*i+1]);
      //   }
      // }
      // else if(m_alpha > 1){
      //   int quot = M/4;
      //   int remainder = M%4;
      //   for (i = 0; i < quot; i++) {
      //     ret[k + 4*i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (4*i)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (4*i) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 4*i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(4*i)+(j+2)]);
      //       }
      //     }
          
      //     ret[k + 4*i+1] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (4*i+1)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (4*i+1) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 4*i+1] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(4*i+1)+(j+2)]);
      //       }
      //     }

      //     ret[k + 4*i + 2] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (4*i+2)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (4*i+2) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 4*i + 2] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(4*i+2)+(j+2)]);
      //       }
      //     }

      //     ret[k + 4*i + 3] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (4*i+3)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (4*i+3) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 4*i + 3] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(4*i+3)+(j+2)]);
      //       }
      //     }
      //   }

      // if(remainder > 0){
      // for(i = 0; i < remainder; i++){
      //   ret[k+4*quot + i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (4*quot+i)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (4*quot+i)+1]);
      //               if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 4*quot + i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(4*quot+i)+(j+2)]);
      //       }
      //     }
      // }
      // }

      // }

//  THE FASTEST!!!
 
      //  if(m_alpha == 1){
      //   for (i = 0; i < M; i++) {
      //     ret[k + i] = SQR(tau) * (SQR(kappa) * fem->doubles[2*i] + fem->doubles[2*i+1]);
      //   }
      // }
      // else if(m_alpha > 1){
      //   int quot = M/3;
      //   int remainder = M%3;
      //   for (i = 0; i < quot; i++) {
      //     ret[k + 3*i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (3*i)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (3*i) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 3*i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(3*i)+(j+2)]);
      //       }
      //     }
          
      //     ret[k + 3*i+1] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (3*i+1)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (3*i+1) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 3*i+1] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(3*i+1)+(j+2)]);
      //       }
      //     }

      //     ret[k + 3*i + 2] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (3*i+2)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (3*i+2) + 1]);

      //     if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 3*i + 2] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(3*i+2)+(j+2)]);
      //       }
      //     }
      //   }
      // if(remainder > 0){
      // for(i = 0; i < remainder; i++){
      //   ret[k+3*quot + i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (3*quot+i)] + m_alpha * 
      //     pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (3*quot+i)+1]);
      //               if(m_alpha>=2){
      //       for(j = 0; j <= (m_alpha-2); j++){
      //         ret[k + 3*quot + i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(3*quot+i)+(j+2)]);
      //       }
      //     }
      // }
      // }

      // } 


    // Testing other splits:

  //  if(m_alpha == 1){
  //       for (i = 0; i < M; i++) {
  //         ret[k + i] = SQR(tau) * (SQR(kappa) * fem->doubles[2*i] + fem->doubles[2*i+1]);
  //       }
  //     }
  //     else if(m_alpha > 1){
  //       int quot = M/6;
  //       int remainder = M%6;
  //       for (i = 0; i < quot; i++) {
  //         ret[k + 6*i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i)+(j+2)]);
  //           }
  //         }
          
  //         ret[k + 6*i+1] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i+1)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i+1) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i+1] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i+1)+(j+2)]);
  //           }
  //         }

  //         ret[k + 6*i + 2] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i+2)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i+2) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i + 2] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i+2)+(j+2)]);
  //           }
  //         }

  //         ret[k + 6*i + 3] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i+3)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i+3) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i + 3] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i+3)+(j+2)]);
  //           }
  //         }

  //         ret[k + 6*i + 4] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i+4)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i+4) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i + 4] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i+4)+(j+2)]);
  //           }
  //          }
          
  //         ret[k + 6*i + 5] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*i+5)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*i+5) + 1]);

  //         if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*i + 5] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*i+5)+(j+2)]);
  //           }
  //         }
  //       }

  //     for(i = 0; i < remainder; i++){
  //       ret[k+6*quot + i] = SQR(tau)*(pow(kappa, 2 * m_alpha) * fem->doubles[(m_alpha+1) * (6*quot+i)] + m_alpha * 
  //         pow(kappa, 2 * (m_alpha-1)) * fem->doubles[(m_alpha+1) * (6*quot+i)+1]);
  //                   if(m_alpha>=2){
  //           for(j = 0; j <= (m_alpha-2); j++){
  //             ret[k + 6*quot + i] += SQR(tau) * (pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2) * fem->doubles[(m_alpha+1)*(6*quot+i)+(j+2)]);
  //           }
  //         }
  //     }

  //     }


      // Fortran implementation
      double sqtau, sqtaukappa, sqtaukappatmp, sqkappatau1, sqkappatau2;
      int one=1;
           sqtau = SQR(tau);
      sqtaukappa = SQR(tau) * SQR(kappa);
      sqkappatau1 = SQR(tau) * SQR(kappa*kappa);
      sqkappatau2 = SQR(tau) * 2 * SQR(kappa);
      dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
      if(m_alpha == 1){
        // double sqtau = SQR(tau);
        // double sqtaukappa = SQR(tau) * SQR(kappa);
        // int one=1;
        // dcopy_(&M, &fem->doubles[M], &one, &ret[k], &one);
        // dscal_(&M, &sqtau, &ret[k], &one); 
        // daxpy_(&M, &sqtaukappa, &fem->doubles[0], &one, &ret[k], &one);

        // int one=1;
        dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
        dscal_(&M, &sqtaukappa, &ret[k], &one); 
        daxpy_(&M, &sqtau, &fem->doubles[M], &one, &ret[k], &one);
      } else if (m_alpha == 2){
        // int one=1;
        dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
        dscal_(&M, &sqkappatau1, &ret[k], &one);
        daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
        daxpy_(&M, &sqtau, &fem->doubles[2*M], &one, &ret[k], &one);
      } else{
        // int one=1;
        dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
        dscal_(&M, &sqkappatau1, &ret[k], &one);
        daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
        if(m_alpha>=2){
          for(j = 2; j<= m_alpha; j++){
            sqkappatau1 = SQR(tau) * pow(kappa, 2 * m_alpha);
            sqkappatau2 = SQR(tau) * m_alpha * pow(kappa, 2 * (m_alpha - 1));
            sqtaukappatmp = SQR(tau) * pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2);
            daxpy_(&M, &sqtaukappatmp, &fem->doubles[j*M], &one, &ret[k], &one);
          }
        }
      }

      // Version using sparsity

      //  double sqtau = SQR(tau);
      // int one=1;
      // if(m_alpha == 1){
      //   // int one=1;
      //   dcopy_(&M, &fem->doubles[0], &one, &ret[k], &one);
      //   dscal_(&M, &sqtau, &ret[k], &one); 
      // } else if (m_alpha == 2){
      //   // int one=1;
      //   dcopy_(&M, &fem->doubles[2*M], &one, &ret[k], &one);
      //   dscal_(&M, &sqtau, &ret[k], &one); 
      //   for(i = 0; i < positions_C->len; i++){
      //     ret[k + positions_G->ints[i]-1] += 2*SQR(tau) * SQR(kappa) * fem->doubles[M+positions_G->ints[i]-1];
      //     ret[k + positions_C->ints[i]-1] += SQR(tau) * SQR(kappa * kappa) * fem->doubles[positions_C->ints[i]-1];
      //   }
      //   for(i = positions_C->len; i<positions_G->len; i++){
      //     ret[k + positions_G->ints[i]-1] += 2*SQR(tau) * SQR(kappa) * fem->doubles[M+positions_G->ints[i]-1];
      //   }

      // } else{
      //   // int one=1;
      // }




    
      // switch(m_alpha){
      //   // double sqtau = SQR(tau);
      //   // double sqtaukappa = SQR(tau) * SQR(kappa);
      //   // int one=1;
      //   // dcopy_(&M, &fem->doubles[M], &one, &ret[k], &one);
      //   // dscal_(&M, &sqtau, &ret[k], &one); 
      //   // daxpy_(&M, &sqtaukappa, &fem->doubles[0], &one, &ret[k], &one);

      //   case 1:
      //   {
      //   dscal_(&M, &sqtaukappa, &ret[k], &one); 
      //   daxpy_(&M, &sqtau, &fem->doubles[M], &one, &ret[k], &one);
      //   break;
      // } case 2:
      // {
      //   dscal_(&M, &sqkappatau1, &ret[k], &one);
      //   daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
      //   daxpy_(&M, &sqtau, &fem->doubles[2*M], &one, &ret[k], &one);
      //   break;
      // } 
      // default:
      // {
      //   dscal_(&M, &sqkappatau1, &ret[k], &one);
      //   daxpy_(&M, &sqkappatau2, &fem->doubles[M], &one, &ret[k], &one);
      //   if(m_alpha>=2){
      //     for(j = 2; j<= m_alpha; j++){
      //       sqkappatau1 = SQR(tau) * pow(kappa, 2 * m_alpha);
      //       sqkappatau2 = SQR(tau) * m_alpha * pow(kappa, 2 * (m_alpha - 1));
      //       sqtaukappatmp = SQR(tau) * pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2);
      //       daxpy_(&M, &sqtaukappatmp, &fem->doubles[j*M], &one, &ret[k], &one);
      //     }
      //   }
      //   break;
      // }
      // }


      // if(m_alpha == 1){
      //   double sqtau = SQR(tau);
      //   double sqtaukappa = SQR(tau) * SQR(kappa);
      //   int one=1;
      //   daxpby_(&M, &sqtaukappa, &fem->doubles[0], &one, &sqtau, &fem->doubles[M], &one, &ret[k]);
      // } else {
      //   int one=1;
      //   double sqkappatau1 = SQR(tau) * pow(kappa, 2 * m_alpha);
      //   double sqkappatau2 = SQR(tau) * m_alpha * pow(kappa, 2 * (m_alpha - 1));
      //   daxpby_(&M, &sqkappatau2, &fem->doubles[M], &one ,&sqkappatau1, &fem->doubles[0], &one, &ret[k]);
      //   if(m_alpha>=2){
      //     for(j = 0; j<= (m_alpha-2); j++){
      //       double sqtaukappatmp = SQR(tau) * pow(kappa, 2*(m_alpha-j-2)) * nChoosek(m_alpha, j+2);
      //       daxpy_(&M, &sqtaukappatmp, &fem->doubles[(j+2)*M], &one, &ret[k], &one);
      //     }
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
      ret[1] = start_ltau;
      ret[2] = start_lkappa;
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

      ret[0] += -0.5 * SQR(lkappa - prior_kappa_meanlog)/(SQR(prior_kappa_sdlog)) - 
      log(prior_kappa_sdlog) - 0.5 * log(2 * M_PI);

      ret[0] += -0.5 * SQR(ltau - prior_tau_meanlog)/(SQR(prior_tau_sdlog)) - 
      log(prior_tau_sdlog) - 0.5 * log(2 * M_PI);

	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}