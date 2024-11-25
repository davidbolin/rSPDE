#include "cgeneric_defs.h"
#include "stdio.h"

// Logit function mapping x in (-1, 1) to (-inf, inf)
double logit_adjusted(double x) {
    if (x <= -1 || x >= 1) {
        fprintf(stderr, "Input x must be in the range (-1, 1).\n");
        return NAN;
    }
    return log((x + 1) / (1 - x));
}

// Inverse logit function mapping (-inf, inf) back to (-1, 1)
double adjusted_inv_logit(double z) {
    return (2.0 / (1.0 + exp(-z))) - 1.0;
}

// Forward transformation: Compute nu from lnu
double forward_nu(double lnu, double nu_upper_bound) {
    return (exp(lnu) / (1.0 + exp(lnu))) * nu_upper_bound;
}

// Inverse transformation: Compute lnu from nu
double inverse_lnu(double nu, double nu_upper_bound) {
    if (nu <= 0 || nu >= nu_upper_bound) {
        fprintf(stderr, "Input nu must be in the range (0, nu_upper_bound).\n");
        return NAN;
    }
    return log(nu / (nu_upper_bound - nu));
}

double *inla_cgeneric_rspde_anisotropic2d_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data) {

    double *ret = NULL;
    int k, i;
    double lhx, lhy, logit_hxy, lnu;
    double lsigma;
    double sigma, hx, hy, hxy;
    char *prior_nu_dist;
    int est_nu = 0;

    // Retrieve parameter values from `data`
    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    //assert(!strcasecmp(data->ints[1]->name, "debug"));
    //int debug = data->ints[1]->ints[0];

    // Basic parameter assertions and retrievals
    assert(!strcasecmp(data->ints[2]->name, "est_nu"));
    est_nu = data->ints[2]->ints[0];

    assert(!strcasecmp(data->ints[3]->name, "rspde_order"));
    int rspde_order = data->ints[3]->ints[0];

    // Retrieve prior means for each parameter
    assert(!strcasecmp(data->doubles[0]->name, "prior.hx.mean"));
    double prior_hx_mean = data->doubles[0]->doubles[0];

    assert(!strcasecmp(data->doubles[1]->name, "prior.hy.mean"));
    double prior_hy_mean = data->doubles[1]->doubles[0];

    assert(!strcasecmp(data->doubles[2]->name, "prior.hxy.mean"));
    double prior_hxy_mean = data->doubles[2]->doubles[0];

    assert(!strcasecmp(data->doubles[3]->name, "prior.sigma.mean"));
    double prior_sigma_mean = data->doubles[3]->doubles[0];

    assert(!strcasecmp(data->doubles[4]->name, "nu"));
    double nu = data->doubles[4]->doubles[0];

    assert(!strcasecmp(data->doubles[5]->name, "start_nu"));
    double start_nu = data->doubles[5]->doubles[0];

    assert(!strcasecmp(data->doubles[6]->name, "prior.nu.loglocation"));
    double prior_nu_loglocation = data->doubles[6]->doubles[0];

    assert(!strcasecmp(data->doubles[7]->name, "prior.nu.mean"));
    double prior_nu_mean = data->doubles[7]->doubles[0];

    assert(!strcasecmp(data->doubles[8]->name, "prior.nu.prec"));
    double prior_nu_prec = data->doubles[8]->doubles[0];

    assert(!strcasecmp(data->doubles[9]->name, "prior.nu.logscale"));
    double prior_nu_logscale = data->doubles[9]->doubles[0];

    assert(!strcasecmp(data->doubles[10]->name, "nu_upper_bound"));
    double nu_upper_bound = data->doubles[10]->doubles[0];

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

    assert(!strcasecmp(data->smats[3]->name, "Hxx"));
    inla_cgeneric_smat_tp *Hxx = data->smats[3];

    assert(!strcasecmp(data->smats[4]->name, "Hyy"));
    inla_cgeneric_smat_tp *Hyy = data->smats[4];

    assert(!strcasecmp(data->smats[5]->name, "Hxy"));
    inla_cgeneric_smat_tp *Hxy = data->smats[5];

    // Retrieve prior precision matrix
    assert(!strcasecmp(data->mats[0]->name, "prior.precision"));
    inla_cgeneric_mat_tp *prior_precision = data->mats[0];

    assert(prior_precision->nrow == 4 && prior_precision->ncol == 4);

    assert(!strcasecmp(data->mats[1]->name, "rational_table"));
    inla_cgeneric_mat_tp *rational_table = data->mats[1];

    if(nu < 0){
        est_nu = 1;
    }

    if (theta) {
        lhx = theta[0];
        lhy = theta[1];
        logit_hxy = theta[2];
        lsigma = theta[3];
        if(est_nu == 1){
            lnu = theta[4];
            nu = forward_nu(lnu, nu_upper_bound);            
        } else {
            lnu = nu = NAN;
        }

        sigma = exp(lsigma);
        hx = exp(lhx);
        hy = exp(lhy);
        hxy = adjusted_inv_logit(logit_hxy);
    } else {   
        lhx = lhy = logit_hxy = lsigma = lnu = sigma = hx = hy = hxy = nu = NAN;
    }

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
            k = 2;
            ret = Calloc(k + M, double);  // Adjust based on Q matrix size
            ret[0] = -1;  // Required value
            ret[1] = M;   
            compute_Q_anisotropic(hx, hy, hxy, sigma, nu, C, Ci, Hxx, Hyy, Hxy, Q, &ret[k], rspde_order, rational_table, est_nu);
            break;

        case INLA_CGENERIC_MU:
        {
            ret = Calloc(1, double);
            ret[0] = 0.0; 
            break;
        }

        case INLA_CGENERIC_INITIAL:
            {
                if(est_nu == 1){
                    ret = Calloc(6, double);
                    ret[0] = 5;  // Number of hyperparameters
                } else{
                    ret = Calloc(5, double);                    
                    ret[0] = 4;
                }

                ret[1] = log(prior_hx_mean);
                ret[2] = log(prior_hy_mean);
                ret[3] = logit_adjusted(prior_hxy_mean);                   
                ret[4] = log(prior_sigma_mean);                   
                if(est_nu == 1){
                    ret[5] = inverse_lnu(start_nu, nu_upper_bound);
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

            double mean_vector[4] = {log(prior_hx_mean), log(prior_hy_mean), logit_adjusted(prior_hxy_mean), log(prior_sigma_mean)};
            double theta_vector[4] = {lhx, lhy, logit_hxy, lsigma};
            ret[0] = logmultnormvdens(4, mean_vector, prior_precision->x, theta_vector);
            
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