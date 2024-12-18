#include "cgeneric_defs.h"
#include "stdio.h"



typedef struct 
{
    double *Q;
    double *mu;
    double *lconst;
    double *theta;
}
my_cache_tp;

double *inla_cgeneric_rspde_fintrinsic_model(inla_cgeneric_cmd_tp cmd, double *theta, 
                                             inla_cgeneric_data_tp *data) {

    double *ret = NULL;
    double *mu_store = NULL;
    double *const_store = NULL;
    double*Q_store = NULL;
    int k, i;
    double lnu, ltau, tau;
    char *prior_nu_dist;
    int est_nu = 0;

    // Retrieve parameter values from `data`
    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    // Basic parameter assertions and retrievals
    assert(!strcasecmp(data->ints[2]->name, "est_nu"));
    est_nu = data->ints[2]->ints[0];

    assert(!strcasecmp(data->ints[3]->name, "rspde_order"));
    int rspde_order = data->ints[3]->ints[0];
    
    assert(!strcasecmp(data->ints[4]->name, "dim"));
    int d = data->ints[4]->ints[0];
    
    assert(!strcasecmp(data->ints[5]->name, "mean_correction"));
    int mean_correction = data->ints[5]->ints[0];

    // Retrieve prior means for each parameter
    assert(!strcasecmp(data->doubles[0]->name, "prior.tau.mean"));
    double prior_tau_mean = data->doubles[0]->doubles[0];
    
    assert(!strcasecmp(data->doubles[1]->name, "prior.tau.precision"));
    double prior_tau_prec = data->doubles[1]->doubles[0];

    assert(!strcasecmp(data->doubles[2]->name, "start_nu"));
    double start_nu = data->doubles[2]->doubles[0];
    
    assert(!strcasecmp(data->doubles[3]->name, "nu"));
    double nu = data->doubles[3]->doubles[0];

    assert(!strcasecmp(data->doubles[4]->name, "prior.nu.loglocation"));
    double prior_nu_loglocation = data->doubles[4]->doubles[0];

    assert(!strcasecmp(data->doubles[5]->name, "prior.nu.mean"));
    double prior_nu_mean = data->doubles[5]->doubles[0];

    assert(!strcasecmp(data->doubles[6]->name, "prior.nu.prec"));
    double prior_nu_prec = data->doubles[6]->doubles[0];

    assert(!strcasecmp(data->doubles[7]->name, "prior.nu.logscale"));
    double prior_nu_logscale = data->doubles[7]->doubles[0];

    assert(!strcasecmp(data->doubles[8]->name, "nu_upper_bound"));
    double nu_upper_bound = data->doubles[8]->doubles[0];
    
    assert(!strcasecmp(data->doubles[9]->name, "scaling"));
    double scaling = data->doubles[9]->doubles[0];

    assert(!strcasecmp(data->chars[2]->name, "prior.nu.dist"));
    prior_nu_dist = &data->chars[2]->chars[0];

    // Retrieve sparse matrix Q
    assert(!strcasecmp(data->smats[0]->name, "Q"));
    inla_cgeneric_smat_tp *Q = data->smats[0];
    int M = Q->n;

    assert(!strcasecmp(data->smats[1]->name, "C"));
    inla_cgeneric_smat_tp *C = data->smats[1];

    assert(!strcasecmp(data->smats[2]->name, "Ci"));
    inla_cgeneric_smat_tp *Ci = data->smats[2];

    assert(!strcasecmp(data->smats[3]->name, "G"));
    inla_cgeneric_smat_tp *G = data->smats[3];
    
    assert(!strcasecmp(data->smats[4]->name, "D"));
    inla_cgeneric_smat_tp *D = data->smats[4];

    assert(!strcasecmp(data->mats[0]->name, "rational_table"));
    inla_cgeneric_mat_tp *rational_table = data->mats[0];

    int n_par = 1;
    if(nu < 0){
        est_nu = 1;
        n_par = 2;
    }

    if (theta) {
        ltau = theta[0];
        tau = exp(ltau);
        if(est_nu == 1){
            lnu = theta[1];
            nu = forward_nu(lnu, nu_upper_bound);            
        } 
    } else {   
        ltau = tau = NAN;
        if(est_nu == 1){
            nu = lnu = NAN;
        }
    }
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
            k = 2;
            ret = Calloc(k + 2 * M, double);
            ret[0] = N;  
            ret[1] = M;  

            for (i = 0; i < M; i++) {
                ret[k++] = Q->i[i];  // Row indices
            }
            for (i = 0; i < M; i++) {
                ret[k++] = Q->j[i];  // Column indices
            }           
            break;
        }

        case INLA_CGENERIC_Q:
        {
            k = 2;
            ret = Calloc(k + M, double);  // Adjust based on Q matrix size
            ret[0] = -1;  // Required value
            ret[1] = M;   
            if(cache) {
                if( (est_nu && cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) ||
                   cache->theta[0] == theta[0]) {
                    memcpy(ret + 2, cache->Q, M * sizeof(double));
                } else {
                    compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, &ret[k], rspde_order, 
                                         rational_table, est_nu, d, 1, 0, 0,
                                         const_store, mu_store, scaling, D);        
                }
            } else {
                compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, &ret[k], rspde_order, 
                                     rational_table, est_nu, d, 1, 0, 0,
                                     const_store, mu_store, scaling, D);    
            }
            break;
        }
            
        case INLA_CGENERIC_MU:
        {
          if (mean_correction) {
            ret = Calloc(1 + N, double);
            ret[0] = N;		/* REQUIRED */
            if(cache) {
                if((est_nu && cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) ||
                   cache->theta[0] == theta[0]) {
                    memcpy(ret + 1, cache->mu, N * sizeof(double));
                } else {
                    Q_store = Calloc(M, double);
                    const_store = Calloc(1, double);
                    compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                         rational_table, est_nu, d, 1, 1, 1,
                                         const_store, &ret[1], scaling, D);  
                    memcpy(cache->Q, Q_store,  M * sizeof(double));
                    memcpy(cache->mu, ret + 1,  N * sizeof(double));
                    memcpy(cache->lconst, const_store,  sizeof(double));
                    memcpy(cache->theta, theta,  n_par * sizeof(double));
                }
            } else {
                Q_store = Calloc(M, double);
                const_store = Calloc(1, double);
                compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                     rational_table, est_nu, d, 1, 1, 1,
                                     const_store, &ret[1], scaling, D);  
                ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                cache->Q = Calloc(M, double);
                cache->mu = Calloc(N, double);
                cache->lconst = Calloc(1, double);
                cache->theta = Calloc(n_par, double);
                memcpy(cache->Q, Q_store,  M * sizeof(double));
                memcpy(cache->mu, ret + 1,  N * sizeof(double));
                memcpy(cache->lconst, const_store,  sizeof(double));
                memcpy(cache->theta, theta,  n_par * sizeof(double));
            }
          } else {
            ret = Calloc(1, double);
            ret[0] = 0.0;
          }
          break;
        }

        case INLA_CGENERIC_INITIAL:
            {
                if(est_nu == 1){
                    ret = Calloc(3, double);
                    ret[0] = 2;  // Number of hyperparameters
                } else{
                    ret = Calloc(2, double);                    
                    ret[0] = 1;
                }
    
                ret[1] = prior_tau_mean;
                if(est_nu == 1){
                    ret[2] = inverse_lnu(start_nu, nu_upper_bound);
                }           
                break;
            }

        case INLA_CGENERIC_LOG_NORM_CONST:
          {
            ret = Calloc(1, double);
            if(cache) {
                if((est_nu && cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) ||
                   cache->theta[0] == theta[0]) {
                    memcpy(ret, cache->lconst, sizeof(double));
                } else {
                    if(mean_correction == 1) {
                        Q_store = Calloc(M, double);
                        mu_store = Calloc(N, double);
                        compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                             rational_table, est_nu, d, 1, 1, 1,
                                             &ret[0], mu_store, scaling, D);  
                        
                        //then cache them and update theta
                        memcpy(cache->Q, Q_store,  M * sizeof(double));
                        memcpy(cache->mu, mu_store,  N * sizeof(double));
                        memcpy(cache->lconst, ret,  sizeof(double));
                        memcpy(cache->theta, theta,  n_par * sizeof(double));
                    } else {
                        Q_store = Calloc(M, double);
                        compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                             rational_table, est_nu, d, 1, 0, 1,
                                             &ret[0], mu_store, scaling, D); 
                        
                        //then cache them and update theta
                        memcpy(cache->Q, Q_store,  M * sizeof(double));
                        memcpy(cache->lconst, ret, sizeof(double));
                        memcpy(cache->theta, theta,  n_par * sizeof(double));
                    } 
                }
            } else {
                if(mean_correction == 1) {
                    Q_store = Calloc(M, double);
                    mu_store = Calloc(N, double);
                    compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                         rational_table, est_nu, d, 1, 1, 1,
                                         &ret[0], mu_store, scaling, D);  
                    
                    //then cache them and update theta
                    ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                    cache->Q = Calloc(M, double);
                    cache->mu = Calloc(N, double);
                    cache->lconst = Calloc(1, double);
                    cache->theta = Calloc(n_par, double);
                    memcpy(cache->Q, Q_store,  M * sizeof(double));
                    memcpy(cache->mu, mu_store,  N * sizeof(double));
                    memcpy(cache->lconst, ret,  sizeof(double));
                    memcpy(cache->theta, theta,  n_par * sizeof(double));
                } else {
                    Q_store = Calloc(M, double);
                    compute_Q_fintrinsic(tau, nu, C, Ci, G, Q, Q_store, rspde_order, 
                                         rational_table, est_nu, d, 1, 0, 1,
                                         &ret[0], mu_store, scaling, D); 
                    
                    //then cache them and update theta
                    ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                    cache->Q = Calloc(M, double);
                    cache->lconst = Calloc(1, double);
                    cache->theta = Calloc(n_par, double);
                    memcpy(cache->Q, Q_store,  M * sizeof(double));
                    memcpy(cache->lconst, ret, sizeof(double));
                    memcpy(cache->theta, theta,  n_par * sizeof(double));
                } 
            }
            break;
          }

        case INLA_CGENERIC_LOG_PRIOR:
        {
            ret = Calloc(1, double);
            
            ret[0] = 0.5 * log(prior_tau_prec) - 0.5 * log(2.0*M_PI);
            ret[0] += 0.5 * SQR(prior_tau_prec * (prior_tau_mean - ltau));

            if(est_nu == 1){
                if(!strcasecmp(prior_nu_dist, "lognormal")){
                  ret[0] += -0.5 * SQR(lnu - prior_nu_loglocation)/(SQR(prior_nu_logscale));
                  ret[0] += -log(prior_nu_logscale) - 0.5 * log(2.0*M_PI);
                  ret[0] -= log(pnorm(log(nu_upper_bound), prior_nu_loglocation, prior_nu_logscale));
                } else { // if(!strcasecmp(prior_nu_dist, "beta")){
                  double s_1 = (prior_nu_mean / nu_upper_bound) * prior_nu_prec;
                  double s_2 = (1 - prior_nu_mean / nu_upper_bound) * prior_nu_prec;
                  ret[0] += logdbeta(nu / nu_upper_bound, s_1, s_2) - log(nu_upper_bound);
                }
            }
            
            break;
        }

        case INLA_CGENERIC_QUIT:
        default:
          break;
        }

    return ret;
}