// cgeneric_aux_hybrid_alpha2.cpp
//
// Eigen-based helpers used by the hybrid Whittle-Matern (alpha = 2)
// cgeneric model defined in cgeneric_rspde_hybrid_alpha2.c.
//
// Two helpers are exposed to C via the cgeneric_defs.h prototypes:
//
//   compute_Q_hybrid_alpha2 : assemble the precision matrix
//                             Q = tau^2 * L * C^{-1} * L,
//                             L = G + kappa^2 * C (mass-lumped C),
//                             and write the values at the graph_i/j
//                             positions (lower-triangular order).
//   compute_mu_hybrid_alpha2 : solve L * mu = X * beta_X with a sparse
//                              Cholesky and return mu.

#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>

extern "C" {
#include "cgeneric_defs.h"
}

#include <vector>

namespace {

// Build the sparse L = G + kappa^2 * C  (with C mass-lumped, supplied
// via its diagonal).
Eigen::SparseMatrix<double> build_L(int N,
                                    const double *C_diag,
                                    const int *G_i, const int *G_j,
                                    const double *G_x, int G_n,
                                    double kappa) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> trips;
    trips.reserve(G_n + N);

    // Add G's triplets.
    for (int k = 0; k < G_n; k++) {
        trips.emplace_back(G_i[k], G_j[k], G_x[k]);
    }
    // Add kappa^2 * C (diagonal).
    const double k2 = kappa * kappa;
    for (int i = 0; i < N; i++) {
        trips.emplace_back(i, i, k2 * C_diag[i]);
    }

    Eigen::SparseMatrix<double> L(N, N);
    L.setFromTriplets(trips.begin(), trips.end());
    return L;
}

// Build C^{-1} as a sparse diagonal (1/C_diag).
Eigen::SparseMatrix<double> build_Cinv(int N, const double *C_diag) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> trips;
    trips.reserve(N);
    for (int i = 0; i < N; i++) {
        double d = C_diag[i];
        if (d <= 0.0) d = 1e-12;
        trips.emplace_back(i, i, 1.0 / d);
    }
    Eigen::SparseMatrix<double> Cinv(N, N);
    Cinv.setFromTriplets(trips.begin(), trips.end());
    return Cinv;
}

}  // namespace

extern "C" void compute_Q_hybrid_alpha2(int N,
                                        const double *C_diag,
                                        const int *G_i, const int *G_j,
                                        const double *G_x, int G_n,
                                        double tau, double kappa,
                                        const int *graph_i,
                                        const int *graph_j,
                                        int M,
                                        double *Q_out) {
    Eigen::SparseMatrix<double> L    = build_L(N, C_diag, G_i, G_j, G_x, G_n, kappa);
    Eigen::SparseMatrix<double> Cinv = build_Cinv(N, C_diag);
    Eigen::SparseMatrix<double> Q    = L * Cinv * L;
    Q = (tau * tau) * Q;

    // Read Q values at the (graph_i, graph_j) positions. Q is
    // mathematically symmetric, so we use .coeff(i, j) which works
    // regardless of which triangle is stored.
    for (int k = 0; k < M; k++) {
        Q_out[k] = Q.coeff(graph_i[k], graph_j[k]);
    }
}

extern "C" void compute_mu_hybrid_alpha2(int N,
                                         const double *C_diag,
                                         const int *G_i, const int *G_j,
                                         const double *G_x, int G_n,
                                         double kappa_mu,
                                         const double *X_x, int p,
                                         const double *beta_x,
                                         double *mu_out) {
    Eigen::SparseMatrix<double> L = build_L(N, C_diag, G_i, G_j, G_x, G_n,
                                            kappa_mu);

    Eigen::VectorXd rhs(N);
    for (int i = 0; i < N; i++) {
        double s = 0.0;
        for (int j = 0; j < p; j++) {
            s += X_x[i * p + j] * beta_x[j];
        }
        rhs[i] = C_diag[i] * s;
    }

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    if (solver.info() != Eigen::Success) {
        for (int i = 0; i < N; i++) mu_out[i] = 0.0;
        return;
    }
    Eigen::VectorXd mu_vec = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
        for (int i = 0; i < N; i++) mu_out[i] = 0.0;
        return;
    }
    for (int i = 0; i < N; i++) mu_out[i] = mu_vec[i];
}
