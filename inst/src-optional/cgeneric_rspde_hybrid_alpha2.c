/* cgeneric_rspde_hybrid_alpha2.c
 *
 * INLA cgeneric implementation of the hybrid Whittle-Matern SPDE model
 * with alpha = 2:
 *
 *     Y(s) = beta_X^T L^{-1} X(s) + tau^{-1} L^{-1} W(s),
 *
 * where L = kappa^2 - Delta. The latent field stored by this cgeneric
 * is the centred process u(s) = tau^{-1} L^{-1} W observed at the FEM
 * mesh nodes, so the precision is
 *
 *     Q = tau^2 * L * C^{-1} * L,
 *
 * and the mu callback returns
 *
 *     mu(s) = beta_X^T (L^{-1} X)(s).
 *
 * Hyperparameters (in the order INLA sees them):
 *
 *   theta[0]              = log(tau)
 *   theta[1]              = log(kappa)
 *   theta[2..(2+p-1)]     = beta_x1, ..., beta_xp   
 */

#include "cgeneric_defs.h"

double *inla_cgeneric_rspde_hybrid_alpha2_model(inla_cgeneric_cmd_tp cmd,
                                                double *theta,
                                                inla_cgeneric_data_tp *data) {

    double *ret = NULL;
    int i, k;

    assert(!strcasecmp(data->ints[0]->name, "n"));
    int N = data->ints[0]->ints[0];
    assert(N > 0);

    assert(!strcasecmp(data->ints[1]->name, "debug"));
    int debug = data->ints[1]->ints[0];
    (void) debug;

    assert(!strcasecmp(data->ints[2]->name, "graph_opt_i"));
    inla_cgeneric_vec_tp *graph_i = data->ints[2];
    int M = graph_i->len;

    assert(!strcasecmp(data->ints[3]->name, "graph_opt_j"));
    inla_cgeneric_vec_tp *graph_j = data->ints[3];
    assert(M == graph_j->len);

    assert(!strcasecmp(data->ints[4]->name, "p"));
    int p = data->ints[4]->ints[0];

    /* separate_kappa_mu: 0 = kappa_mu tied to kappa (one fewer hyperpar),
     * 1 = kappa_mu is its own hyperparameter slotted between
     * log(kappa) and beta_x1.
     */
    assert(!strcasecmp(data->ints[5]->name, "separate_kappa_mu"));
    int separate_kappa_mu = data->ints[5]->ints[0];

    assert(!strcasecmp(data->doubles[0]->name, "C_diag"));
    inla_cgeneric_vec_tp *C_diag_v = data->doubles[0];
    assert(C_diag_v->len == N);

    assert(!strcasecmp(data->smats[0]->name, "G_mat"));
    inla_cgeneric_smat_tp *G_in = data->smats[0];

    assert(!strcasecmp(data->mats[0]->name, "X_mat"));
    inla_cgeneric_mat_tp *X_in = data->mats[0];
    assert(X_in->nrow == N);
    assert(X_in->ncol == p);

    assert(!strcasecmp(data->doubles[1]->name, "theta.prior.mean"));
    inla_cgeneric_vec_tp *theta_prior_mean = data->doubles[1];

    assert(!strcasecmp(data->mats[1]->name, "theta.prior.prec"));
    inla_cgeneric_mat_tp *theta_prior_prec = data->mats[1];

    assert(!strcasecmp(data->doubles[2]->name, "start.theta"));
    inla_cgeneric_vec_tp *start_theta = data->doubles[2];

    assert(!strcasecmp(data->doubles[3]->name, "beta_x.prior.mean"));
    inla_cgeneric_vec_tp *beta_x_prior_mean = data->doubles[3];

    assert(!strcasecmp(data->doubles[4]->name, "beta_x.prior.prec"));
    inla_cgeneric_vec_tp *beta_x_prior_prec = data->doubles[4];

    /* Optional kappa_mu prior (only used when separate_kappa_mu = 1). */
    assert(!strcasecmp(data->doubles[5]->name, "kappa_mu.prior.mean"));
    inla_cgeneric_vec_tp *kappa_mu_prior_mean = data->doubles[5];

    assert(!strcasecmp(data->doubles[6]->name, "kappa_mu.prior.prec"));
    inla_cgeneric_vec_tp *kappa_mu_prior_prec = data->doubles[6];

    /* Decode theta:
     *   theta[0]            = log(tau)
     *   theta[1]            = log(kappa)
     *   [if separate]
     *     theta[2]          = log(kappa_mu)
     *     theta[3..]        = beta_x
     *   [else]
     *     theta[2..]        = beta_x
     */
    double ltau, lkappa, lkappa_mu, tau, kappa, kappa_mu;
    double *beta_x = (double*) Calloc(p > 0 ? p : 1, double);
    if (theta) {
        ltau   = theta[0];
        lkappa = theta[1];
        tau    = exp(ltau);
        kappa  = exp(lkappa);
        if (separate_kappa_mu) {
            lkappa_mu = theta[2];
            kappa_mu  = exp(lkappa_mu);
            for (i = 0; i < p; i++) beta_x[i] = theta[3 + i];
        } else {
            lkappa_mu = lkappa;
            kappa_mu  = kappa;
            for (i = 0; i < p; i++) beta_x[i] = theta[2 + i];
        }
    } else {
        ltau = lkappa = lkappa_mu = NAN;
        tau = kappa = kappa_mu = NAN;
        for (i = 0; i < p; i++) beta_x[i] = NAN;
    }

    switch (cmd) {

    case INLA_CGENERIC_VOID:
        assert(!(cmd == INLA_CGENERIC_VOID));
        break;

    case INLA_CGENERIC_GRAPH:
        {
            k = 2;
            ret = Calloc(k + 2 * M, double);
            ret[0] = N;
            ret[1] = M;
            for (i = 0; i < M; i++) ret[k++] = graph_i->ints[i];
            for (i = 0; i < M; i++) ret[k++] = graph_j->ints[i];
            break;
        }

    case INLA_CGENERIC_Q:
        {
            k = 2;
            ret = Calloc(k + M, double);
            ret[0] = -1;        /* REQUIRED */
            ret[1] = M;
            compute_Q_hybrid_alpha2(N,
                                    C_diag_v->doubles,
                                    G_in->i, G_in->j, G_in->x, G_in->n,
                                    tau, kappa,
                                    graph_i->ints, graph_j->ints, M,
                                    &ret[k]);
            break;
        }

    case INLA_CGENERIC_MU:
        {
            /* If all beta_x are exactly zero, return the "no mean"
             * signal so that INLA falls back to its default zero-mean
             * handling: an explicit [N, 0, ..., 0] return is not
             * numerically equivalent in INLA's marginal-likelihood
             * code, and would create a spurious normalisation
             * difference relative to a standard (zero-mean) latent
             * model.
             */
            int all_zero = 1;
            for (i = 0; i < p; i++) {
                if (beta_x[i] != 0.0) { all_zero = 0; break; }
            }
            if (all_zero) {
                ret = Calloc(1, double);
                ret[0] = 0.0;
                break;
            }
            ret = Calloc(1 + N, double);
            ret[0] = N;         /* REQUIRED */
            compute_mu_hybrid_alpha2(N,
                                     C_diag_v->doubles,
                                     G_in->i, G_in->j, G_in->x, G_in->n,
                                     kappa_mu,
                                     X_in->x, p, beta_x,
                                     &ret[1]);
            break;
        }

    case INLA_CGENERIC_INITIAL:
        {
            int n_par = 2 + (separate_kappa_mu ? 1 : 0) + p;
            ret = Calloc(1 + n_par, double);
            ret[0] = n_par;
            ret[1] = start_theta->doubles[0];   /* log(tau)   */
            ret[2] = start_theta->doubles[1];   /* log(kappa) */
            int off = 2;
            if (separate_kappa_mu) {
                ret[3] = start_theta->doubles[2];   /* log(kappa_mu) */
                off = 3;
            }
            for (i = 0; i < p; i++) {
                ret[1 + off + i] = start_theta->doubles[off + i];
            }
            break;
        }

    case INLA_CGENERIC_LOG_NORM_CONST:
        /* let INLA compute it from Q */
        break;

    case INLA_CGENERIC_LOG_PRIOR:
        {
            ret = Calloc(1, double);
            ret[0] = 0.0;

            /* (log(tau), log(kappa)) ~ N(theta_prior_mean, theta_prior_prec^{-1}) */
            double th2[2];
            th2[0] = ltau;
            th2[1] = lkappa;
            ret[0] += logmultnormvdens(2, theta_prior_mean->doubles,
                                       theta_prior_prec->x, th2);

            /* log(kappa_mu) ~ N(mean, 1/prec) when separate */
            if (separate_kappa_mu) {
                double m_km = kappa_mu_prior_mean->doubles[0];
                double pr_km = kappa_mu_prior_prec->doubles[0];
                if (pr_km <= 0.0) pr_km = 1e-6;
                double dx_km = lkappa_mu - m_km;
                ret[0] += 0.5 * log(pr_km) - 0.5 * log(2.0 * M_PI)
                          - 0.5 * pr_km * dx_km * dx_km;
            }

            /* beta_x_j ~ N(mean_j, 1/prec_j), independent */
            for (i = 0; i < p; i++) {
                double m_j = beta_x_prior_mean->doubles[i];
                double pr  = beta_x_prior_prec->doubles[i];
                if (pr <= 0.0) pr = 1e-6;
                double dx = beta_x[i] - m_j;
                ret[0] += 0.5 * log(pr) - 0.5 * log(2.0 * M_PI)
                          - 0.5 * pr * dx * dx;
            }
            break;
        }

    case INLA_CGENERIC_QUIT:
    default:
        break;
    }

    free(beta_x);
    return ret;
}
