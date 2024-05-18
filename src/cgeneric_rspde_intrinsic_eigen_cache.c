#include "cgeneric_defs.h"
#include "cgeneric.h"
#include "stdio.h"
//#include "omp.h"
    // #include "gsl/gsl_vector_double.h"
    
    typedef struct 
{
    double *Q;
    double *mu;
    double *lconst;
    double *theta;
}
my_cache_tp;

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_rspde_intrinsic_eigen_cache(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {
    
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
                        
                if(cache) {
                    // cache exists, check if theta is the same
                    if(cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) {
                        //printf("Q: reusing Q for cache_idx = %1d\n", cache_idx);
                        // same parameters, return Q
                        memcpy(ret + 2, cache->Q, M * sizeof(double));
                    } else {
                        //printf("Q: Computing Q for cache_idx = %1d\n", cache_idx);
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
                    }
                } else {
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
                }
                break;
            }
        
        case INLA_CGENERIC_MU:
            {
                if(mean_correction == 1) { 
                    ret = Calloc(1 + n_mesh, double);
                    ret[0] = n_mesh;		/* REQUIRED */
                        
                        
                        if(cache) {
                            // cache exists, check if theta is the same
                            if(cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) {
                                // same parameters, return mu
                               //printf("mu: reusing mu for cache_idx = %1d\n", cache_idx);
                                memcpy(ret + 1, cache->mu, n_mesh * sizeof(double));
                            } else {
                                // new parameters, compute the quantities
                                //printf("mu: Computing mu for cache_idx = %1d\n", cache_idx);
                                if(true_scaling == 1) {
                                    Q_store = Calloc(M, double);
                                    const_store = Calloc(1, double);
                                    compute_Q_intrinsic(n_mesh, 
                                                        C->x, C->i, C->j, C->n,
                                                        G->x, G->i, G->j, G->n,
                                                        theta, 
                                                        Q_store,
                                                        alpha,
                                                        1,
                                                        1,
                                                        1,
                                                        const_store, 
                                                        &ret[1]);    
                                    //then cache them and update theta
                                    memcpy(cache->Q, Q_store,  M * sizeof(double));
                                    memcpy(cache->mu, ret + 1,  n_mesh * sizeof(double));
                                    memcpy(cache->lconst, const_store,  sizeof(double));
                                    memcpy(cache->theta, theta,  2 * sizeof(double));
                                } else {
                                    Q_store = Calloc(M, double);
                                    compute_Q_intrinsic(n_mesh, 
                                                        C->x, C->i, C->j, C->n,
                                                        G->x, G->i, G->j, G->n,
                                                        theta, 
                                                        Q_store,
                                                        alpha,
                                                        1,
                                                        1,
                                                        0,
                                                        const_store, 
                                                        &ret[1]);    
                                    //then cache them and update theta
                                    memcpy(cache->Q, Q_store,  M * sizeof(double));
                                    memcpy(cache->mu, ret + 1,  n_mesh * sizeof(double));
                                    memcpy(cache->theta, theta,  2 * sizeof(double));
                                } 
                            }
                        } else {
                            // cache empty, compute quantities
                            if(true_scaling == 1) {
                                Q_store = Calloc(M, double);
                                const_store = Calloc(1, double);
                                compute_Q_intrinsic(n_mesh, 
                                                    C->x, C->i, C->j, C->n,
                                                    G->x, G->i, G->j, G->n,
                                                    theta, 
                                                    Q_store,
                                                    alpha,
                                                    1,
                                                    1,
                                                    1,
                                                    const_store, 
                                                    &ret[1]);    
                                //then allocate memory in the cache and store it and theta
                                ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                                cache->Q = Calloc(M, double);
                                cache->mu = Calloc(n_mesh, double);
                                cache->lconst = Calloc(1, double);
                                cache->theta = Calloc(2, double);
                                memcpy(cache->Q, Q_store,  M * sizeof(double));
                                memcpy(cache->mu, ret + 1,  n_mesh * sizeof(double));
                                memcpy(cache->lconst, const_store,  sizeof(double));
                                memcpy(cache->theta, theta,  2 * sizeof(double));
                            } else  {
                                Q_store = Calloc(M, double);
                                compute_Q_intrinsic(n_mesh, 
                                                    C->x, C->i, C->j, C->n,
                                                    G->x, G->i, G->j, G->n,
                                                    theta, 
                                                    Q_store,
                                                    alpha,
                                                    1,
                                                    1,
                                                    0,
                                                    const_store, 
                                                    &ret[1]);   
                                //then allocate memory in the cache and store it and theta
                                ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                                cache->Q = Calloc(M, double);
                                cache->mu = Calloc(n_mesh, double);
                                cache->theta = Calloc(2, double);
                                memcpy(cache->Q, Q_store,  M * sizeof(double));
                                memcpy(cache->mu, ret + 1,  n_mesh * sizeof(double));
                                memcpy(cache->theta, theta,  2 * sizeof(double));
                            } 
                        }
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
                    if(cache) {
                        // cache exists, check if theta is the same
                        if(cache->theta[0] == theta[0] && cache->theta[1] == theta[1]) {
                            // same parameters, return const
                            //printf("const: reusing const for cache_idx = %1d\n", cache_idx);
                            memcpy(ret, cache->lconst, sizeof(double));
                        } else {
                            // new parameters, compute the quantities
                            //printf("const: Computing const for cache_idx = %1d\n", cache_idx);
                            if(mean_correction == 1) {
                                Q_store = Calloc(M, double);
                                mu_store = Calloc(n_mesh, double);
                                compute_Q_intrinsic(n_mesh, 
                                                    C->x, C->i, C->j, C->n,
                                                    G->x, G->i, G->j, G->n,
                                                    theta, 
                                                    Q_store,
                                                    alpha,
                                                    1,
                                                    1,
                                                    1,
                                                    &ret[0], 
                                                    mu_store);    
                                //then cache them and update theta
                                memcpy(cache->Q, Q_store,  M * sizeof(double));
                                memcpy(cache->mu, mu_store,  n_mesh * sizeof(double));
                                memcpy(cache->lconst, ret,  sizeof(double));
                                memcpy(cache->theta, theta,  2 * sizeof(double));
                            } else {
                                Q_store = Calloc(M, double);
                                compute_Q_intrinsic(n_mesh, 
                                                    C->x, C->i, C->j, C->n,
                                                    G->x, G->i, G->j, G->n,
                                                    theta, 
                                                    Q_store,
                                                    alpha,
                                                    1,
                                                    0,
                                                    1,
                                                    &ret[0], 
                                                    mu_store);    
                                //then cache them and update theta
                                memcpy(cache->Q, Q_store,  M * sizeof(double));
                                memcpy(cache->lconst, ret, sizeof(double));
                                memcpy(cache->theta, theta,  2 * sizeof(double));
                            } 
                        }
                    } else {
                        // cache empty, compute quantities
                        if(mean_correction == 1) {
                            Q_store = Calloc(M, double);
                            mu_store = Calloc(n_mesh, double);
                            compute_Q_intrinsic(n_mesh, 
                                                C->x, C->i, C->j, C->n,
                                                G->x, G->i, G->j, G->n,
                                                theta, 
                                                Q_store,
                                                alpha,
                                                1,
                                                1,
                                                1,
                                                &ret[0], 
                                                mu_store);    
                            //then allocate memory in the cache and store it and theta
                            ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                            cache->Q = Calloc(M, double);
                            cache->mu = Calloc(n_mesh, double);
                            cache->lconst = Calloc(1, double);
                            cache->theta = Calloc(2, double);
                            memcpy(cache->Q, Q_store,  M * sizeof(double));
                            memcpy(cache->mu, mu_store,  n_mesh * sizeof(double));
                            memcpy(cache->lconst, ret,  sizeof(double));
                            memcpy(cache->theta, theta,  2 * sizeof(double));
                        } else  {
                            Q_store = Calloc(M, double);
                            compute_Q_intrinsic(n_mesh, 
                                                C->x, C->i, C->j, C->n,
                                                G->x, G->i, G->j, G->n,
                                                theta, 
                                                Q_store,
                                                alpha,
                                                1,
                                                0,
                                                1,
                                                &ret[0], 
                                                mu_store);   
                            //then allocate memory in the cache and store it and theta
                            ((my_cache_tp **) data->cache)[cache_idx] = cache = Calloc(1, my_cache_tp);
                            cache->Q = Calloc(M, double);
                            cache->lconst = Calloc(1, double);
                            cache->theta = Calloc(2, double);
                            memcpy(cache->Q, Q_store,  M * sizeof(double));
                            memcpy(cache->lconst, ret,  sizeof(double));
                            memcpy(cache->theta, theta,  2 * sizeof(double));
                        } 
                    }
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