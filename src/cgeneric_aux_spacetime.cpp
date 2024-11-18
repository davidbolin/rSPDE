#include <Eigen/Sparse>
#include "cgeneric_defs.h"
#include "cgeneric_cpp.h"

extern "C" {

double nChoosek(int n, int k);

void compute_Q_dim1(
    double kappa, double sigma, double gamma, double rho, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list
);
void compute_Q_dim2(
    double kappa, double sigma, double gamma, double rho_1, double rho_2, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list
);

} // extern "C"

// Define a type alias for SparseMatrix with ColMajor storage
typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrixColMajor;

// Helper function to convert inla_cgeneric_smat_tp to Eigen::SparseMatrix with ColMajor storage
SparseMatrixColMajor convertInlaToEigen(const inla_cgeneric_smat_tp* inlaMat) {
    if (!inlaMat || inlaMat->n <= 0) {
        throw std::invalid_argument("Invalid INLA sparse matrix");
    }
    SparseMatrixColMajor eigenMat(inlaMat->nrow, inlaMat->ncol);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(inlaMat->n);

    for (int idx = 0; idx < inlaMat->n; ++idx) {
        triplets.emplace_back(inlaMat->i[idx], inlaMat->j[idx], inlaMat->x[idx]);
    }

    eigenMat.setFromTriplets(triplets.begin(), triplets.end());
    return eigenMat;
}

// Function to create matrix L based on Glist, kappa, and n, using ColMajor storage
SparseMatrixColMajor make_L(int n, double kappa, inla_cgeneric_smat_tp **Glist) {
    SparseMatrixColMajor L = convertInlaToEigen(Glist[0]);
    L.setZero();  // Start with L as a zero matrix of the same size

    for (int k = 0; k <= n; k++) {
        double coeff = nChoosek(n, k) * pow(kappa, 2 * (n - k));
        SparseMatrixColMajor Gk = convertInlaToEigen(Glist[k]);

        L += coeff * Gk;
    }

    return L;
}

// Main function to compute Q for the 1D case (d = 1)
void compute_Q_dim1(
    double kappa, double sigma, double gamma, double rho, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist,
    inla_cgeneric_smat_tp** Ctlist,
    inla_cgeneric_smat_tp** B0list,
    inla_cgeneric_smat_tp*** M2list
) {
    // Initialize Q with the first term using ColMajor storage
    SparseMatrixColMajor Q = make_L(beta, kappa, Gtlist) + 2 * gamma * make_L(beta + alpha, kappa, B0list);

    // Loop over alpha to add the contributions from each term
    for (int k = 0; k <= alpha; ++k) {
        double gamma_sq_rho = gamma * gamma * nChoosek(alpha, k) * pow(rho, 2 * k);

        // Compute the gamma^2 term
        SparseMatrixColMajor Ct_sum = make_L(beta + 2 * (alpha - k), kappa, &Ctlist[k]);

        // Add this term to Q
        Q += gamma_sq_rho * Ct_sum;

        // Calculate the M2 term
        SparseMatrixColMajor M2_term = make_L(beta + alpha - k, kappa, M2list[k]);

        // Add the M2 term to Q
        double factor = -0.5 * pow(-1, k / 2) * gamma * nChoosek(alpha, k) * (1 - pow(-1, k)) * pow(rho, k);

        // Ensure that M2_term and its transpose have the same storage order
        SparseMatrixColMajor M2_term_transpose = M2_term.transpose();

        Q += factor * (M2_term + M2_term_transpose);
    }

    Q /= sigma * sigma;

    // Extract the values from Q into the result array, using only lower triangular part
    Eigen::SparseMatrix<double, Eigen::ColMajor> Q_triang = Q.triangularView<Eigen::Lower>();

    int count = 0;
    for (int k = 0; k < Q_triang.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(Q_triang, k); it; ++it) {
            result[count] = it.value();
            count++;
        }
    }
}

// Main function to compute Q for the 2D case (d = 2)
void compute_Q_dim2(
    double kappa, double sigma, double gamma, double rho_1, double rho_2, int beta, int alpha,
    double* result,  // Output array for Q
    inla_cgeneric_smat_tp** Gtlist, 
    inla_cgeneric_smat_tp** Ctlist, 
    inla_cgeneric_smat_tp** B0list, 
    inla_cgeneric_smat_tp*** M2list
) {
    SparseMatrixColMajor Q = make_L(beta, kappa, Gtlist) + 2 * gamma * make_L(beta + alpha, kappa, B0list);

    // Add the main Ctlist term
    Q += gamma * gamma * make_L(beta + 2 * alpha, kappa, Ctlist);

    if (alpha == 1) {
        SparseMatrixColMajor tmp = rho_1 * rho_1 * make_L(beta, kappa, M2list[0]) +
                                   rho_2 * rho_2 * make_L(beta, kappa, M2list[1]) +
                                   2 * rho_1 * rho_2 * make_L(beta, kappa, M2list[2]);

        SparseMatrixColMajor tmp_transp = tmp.transpose();
        
        Q += 0.5 * gamma * gamma * (tmp + tmp_transp);

        // Calculate M2 term
        SparseMatrixColMajor M2 = rho_1 * make_L(beta, kappa, M2list[3]) +
                                  rho_2 * make_L(beta, kappa, M2list[4]);

        SparseMatrixColMajor M2_transp = M2.transpose();
        
        Q -= gamma * (M2 + M2_transp);
    }

    // Final scaling by sigma squared
    Q /= sigma * sigma;

    // Extract the values from Q into the result array
    Eigen::SparseMatrix<double> Q_triang(Q.rows(), Q.rows());
    Q_triang = Q.triangularView<Eigen::Lower>();    

    int count = 0;
    for (int k = 0; k < Q_triang.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang, k); it; ++it) {
            result[count] = it.value();
            count++;
        }
    }
}