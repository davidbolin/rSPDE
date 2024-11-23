#include "cgeneric_defs.h"
#include "stdio.h"

double *inla_cgeneric_rspde_spacetime_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data) {

    double *ret = NULL;
    int k, i;
    double lkappa, lsigma, lgamma, rho, rho2;
    double kappa, sigma, gamma;
    double prior_rho2_mean;

    // Retrieve parameter values from `data`
    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    //assert(!strcasecmp(data->ints[1]->name, "debug"));
    //int debug = data->ints[1]->ints[0];

    // Basic parameter assertions and retrievals
    assert(!strcasecmp(data->ints[2]->name, "d"));
    int d = data->ints[2]->ints[0];
    
    assert(!strcasecmp(data->ints[3]->name, "n_Gtlist"));
    int n_Gtlist = data->ints[3]->ints[0];
    
    assert(!strcasecmp(data->ints[4]->name, "n_Ctlist"));
    int n_Ctlist = data->ints[4]->ints[0];
    
    assert(!strcasecmp(data->ints[5]->name, "n_B0list"));
    int n_B0list = data->ints[5]->ints[0];
    
    assert(!strcasecmp(data->ints[6]->name, "n_M2list"));
    int n_M2list = data->ints[6]->ints[0];
    
    //assert(!strcasecmp(data->ints[7]->name, "n_M2list2"));
    //int n_M2list2 = data->ints[7]->ints[0];
    
    assert(!strcasecmp(data->ints[8]->name, "drift"));
    int drift = data->ints[8]->ints[0];

    // Retrieve prior means for each parameter
    assert(!strcasecmp(data->doubles[0]->name, "prior.kappa.mean"));
    double prior_kappa_mean = data->doubles[0]->doubles[0];

    assert(!strcasecmp(data->doubles[1]->name, "prior.sigma.mean"));
    double prior_sigma_mean = data->doubles[1]->doubles[0];

    assert(!strcasecmp(data->doubles[2]->name, "prior.gamma.mean"));
    double prior_gamma_mean = data->doubles[2]->doubles[0];

    assert(!strcasecmp(data->doubles[3]->name, "prior.rho.mean"));
    double prior_rho_mean = data->doubles[3]->doubles[0];

    if(d == 2){
        prior_rho2_mean = data->doubles[3]->doubles[1];
    } else {
        prior_rho2_mean = NAN;
    }

    // Beta and Alpha
    assert(!strcasecmp(data->doubles[4]->name, "beta"));
    double beta = data->doubles[4]->doubles[0];

    assert(!strcasecmp(data->doubles[5]->name, "alpha"));
    double alpha = data->doubles[5]->doubles[0];

    // Retrieve sparse matrix Q
    assert(!strcasecmp(data->smats[0]->name, "Q"));
    inla_cgeneric_smat_tp *Q = data->smats[0];
    int M = Q->n;

    int smat_index = 1;

    // Retrieve prior precision matrix
    assert(!strcasecmp(data->mats[0]->name, "prior.precision"));
    inla_cgeneric_mat_tp *prior_precision = data->mats[0];

    if (d == 2) {
        if(drift == 1){
            assert(prior_precision->nrow == 5 && prior_precision->ncol == 5);
        } else{
            assert(prior_precision->nrow == 4 && prior_precision->ncol == 4);
        }
    } else if (d == 1) {
        if(drift == 1){
            assert(prior_precision->nrow == 4 && prior_precision->ncol == 4);
        } else{
            assert(prior_precision->nrow == 3 && prior_precision->ncol == 3);
        }
    } else {
        // Handle unexpected dimension cases if necessary
        assert(0 && "Unsupported dimension for prior precision matrix");
    }

    // Process list of matrices (Gtlist, Ctlist, B0list, M2list, M2list2)
    inla_cgeneric_smat_tp **Gtlist = malloc(n_Gtlist * sizeof(inla_cgeneric_smat_tp *));
    for (int i = 0; i < n_Gtlist; i++) {
        char expected_name[20];
        sprintf(expected_name, "Gtlist%d", i + 1);
        assert(!strcasecmp(data->smats[smat_index]->name, expected_name));
        Gtlist[i] = data->smats[smat_index++];
    }

    inla_cgeneric_smat_tp **Ctlist = malloc(n_Ctlist * sizeof(inla_cgeneric_smat_tp *));
    for (int i = 0; i < n_Ctlist; i++) {
        char expected_name[20];
        sprintf(expected_name, "Ctlist%d", i + 1);
        assert(!strcasecmp(data->smats[smat_index]->name, expected_name));
        Ctlist[i] = data->smats[smat_index++];
    }

    inla_cgeneric_smat_tp **B0list = malloc(n_B0list * sizeof(inla_cgeneric_smat_tp *));
    for (int i = 0; i < n_B0list; i++) {
        char expected_name[20];
        sprintf(expected_name, "B0list%d", i + 1);
        assert(!strcasecmp(data->smats[smat_index]->name, expected_name));
        B0list[i] = data->smats[smat_index++];
    }

    // Two-level structure for M2list and M2list2
    inla_cgeneric_smat_tp ***M2list = malloc(n_M2list * sizeof(inla_cgeneric_smat_tp **));
    for (int i = 0; i < n_M2list; i++) {
        char length_name[20];
        sprintf(length_name, "n_M2list_%d", i + 1);
        
        int n_M2list_i = 0;
        for (int j = 0; j < data->n_ints; j++) {
            if (!strcasecmp(data->ints[j]->name, length_name)) {
                n_M2list_i = data->ints[j]->ints[0];
                break;
            }
        }
        assert(n_M2list_i > 0);

        M2list[i] = malloc(n_M2list_i * sizeof(inla_cgeneric_smat_tp *));
        for (int j = 0; j < n_M2list_i; j++) {
            char element_name[30];
            sprintf(element_name, "M2list%d_%d", i + 1, j + 1);
            assert(!strcasecmp(data->smats[smat_index]->name, element_name));
            M2list[i][j] = data->smats[smat_index++];
        }
    }

    if (theta) {
        lkappa = theta[0];
        lsigma = theta[1];
        lgamma = theta[2];
        if(drift == 1){
            rho = theta[3];  // No exponential transformation for rho
            if(d == 2){
                rho2 = theta[4];
            } else{
                rho2 = 0.0;
            }
        } else{
            rho = 0.0;
            rho2 = 0.0;
        }

        kappa = exp(lkappa);
        sigma = exp(lsigma);
        gamma = exp(lgamma);
    } else {   
        lkappa = lsigma = lgamma = kappa = sigma = gamma = rho = rho2 = NAN;
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

            if(d == 1){
                compute_Q_dim1(kappa, sigma, gamma, rho, beta, alpha, &ret[k], Gtlist, 
                        Ctlist, B0list, M2list);
            } else if(d == 2){
                compute_Q_dim2(kappa, sigma, gamma, rho, rho2, beta, alpha, &ret[k], Gtlist, 
                        Ctlist, B0list, M2list);
            }

            break;

        case INLA_CGENERIC_MU:
        {
            ret = Calloc(1, double);
            ret[0] = 0.0; 
            break;
        }

        case INLA_CGENERIC_INITIAL:
            {
                if(drift == 1){
                    if(d == 1){
                        ret = Calloc(5, double);
                        ret[0] = 4;  // Number of hyperparameters
                    } else{
                        ret = Calloc(6, double);
                        ret[0] = 5;
                    }
                } else{
                    ret = Calloc(4, double);                    
                    ret[0] = 3;
                }

                ret[1] = log(prior_kappa_mean);                   
                ret[2] = log(prior_sigma_mean);                   
                ret[3] = log(prior_gamma_mean);        
                if(drift == 1){
                    ret[4] = prior_rho_mean;                     
                    if(d == 2){
                        ret[5] = prior_rho2_mean;
                    } 
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
            
            if(drift == 1){
                if(d == 1){
                    double mean_vector[4] = {log(prior_kappa_mean), log(prior_sigma_mean), log(prior_gamma_mean), prior_rho_mean};
                    double theta_vector[4] = {lkappa, lsigma, lgamma, rho};
                    ret[0] = logmultnormvdens(4, mean_vector, prior_precision->x, theta_vector);
                } else{
                    double mean_vector[5] = {log(prior_kappa_mean), log(prior_sigma_mean), log(prior_gamma_mean), prior_rho_mean, prior_rho2_mean};
                    double theta_vector[5] = {lkappa, lsigma, lgamma, rho, rho2};
                    ret[0] = logmultnormvdens(5, mean_vector, prior_precision->x, theta_vector);
                }
            } else{
                double mean_vector[3] = {log(prior_kappa_mean), log(prior_sigma_mean), log(prior_gamma_mean)};
                double theta_vector[3] = {lkappa, lsigma, lgamma};
                ret[0] = logmultnormvdens(3, mean_vector, prior_precision->x, theta_vector);
            }
            
            break;
        }

        case INLA_CGENERIC_QUIT:
        default:
          break;
        }

    return ret;
}