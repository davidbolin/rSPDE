#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include "cgeneric_defs.h"
#include "cgeneric_cpp.h"

/**
 * Convert an inla_cgeneric_mat_tp to Eigen::MatrixXd.
 * @param mat The input matrix in the inla_cgeneric_mat_tp format.
 * @return The converted Eigen::MatrixXd.
 */
Eigen::MatrixXd inlaToEigenMatrix(const inla_cgeneric_mat_tp *mat) {
    if (mat == nullptr || mat->x == nullptr) {
        throw std::invalid_argument("Input matrix or data pointer is null.");
    }

    Eigen::MatrixXd eigenMat(mat->nrow, mat->ncol);

    // Map row-major data to column-major Eigen::MatrixXd
    for (int i = 0; i < mat->nrow; ++i) {
        for (int j = 0; j < mat->ncol; ++j) {
            eigenMat(i, j) = mat->x[i * mat->ncol + j];
        }
    }

    return eigenMat;
}

Eigen::SparseMatrix<double> anisotropic_precision(
    const Eigen::SparseMatrix<double>& L,  // Input matrix L (FEM operator)
    double tau,                            // Scaling factor (tau)
    const Eigen::SparseMatrix<double>& C,  // Mass matrix C
    const Eigen::SparseMatrix<double>& Ci,  // Inverse Mass matrix C
    const Eigen::SparseMatrix<double>& CiL,  // Ci * L
    double alpha,                          // Fractional exponent alpha
    int m_alpha,                           // Integer part of alpha
    int rspde_order,                       // Order of rational approximation
    double scale_factor,                   // Scaling factor for precision
    const Eigen::MatrixXd& rational_table
) {
    int size = L.rows(); // Size of the individual blocks
    int total_size = (rspde_order + 1) * size;

    Eigen::SparseMatrix<double> Q(total_size, total_size); // Resulting matrix
    Eigen::SparseMatrix<double> Q_tmp(size, size); // Temporary matrix for each block

    int row_nu = static_cast<int>(std::round(1000 * cut_decimals(alpha)));
    Eigen::VectorXd rat_r = rational_table.row(row_nu).segment(1, rspde_order);
    Eigen::VectorXd rat_p = rational_table.row(row_nu).segment(1 + rspde_order, rspde_order);
    double k_rat = rational_table(row_nu, 1 + 2 * rspde_order);

    Eigen::SparseMatrix<double> M = Eigen::SparseMatrix<double>(L.rows(), L.cols());
    M.setIdentity();    

    // Add the integer part contribution to M
    if (m_alpha > 0) {
        for (int i = 0; i < m_alpha; ++i) {
            M = M * CiL;
        }
    }

    // Compute the fractional part of the precision matrix
    for (int k = 0; k < rspde_order; ++k) {
        Eigen::SparseMatrix<double> Q_tmp = (L - rat_p(k) * C) * M / rat_r(k);
        
        // Insert Q_tmp into the appropriate block of Q
        for (int m = 0; m < Q_tmp.outerSize(); ++m) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Q_tmp, m); it; ++it) {
                Q.insert(it.row() + size * k, it.col() + size * k) = it.value();
            }
        }
    }

    // Add the k-part to the matrix
    Eigen::SparseMatrix<double> K;
    if (m_alpha > 0) {
        K = C;
        for (int i = 0; i < m_alpha; ++i) {
            K = K * CiL;
        }
    } else {
        K = Ci;
    }
    K = K / k_rat;

    for (int m = 0; m < K.outerSize(); ++m) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, m); it; ++it) {
            Q.insert(it.row() + size * rspde_order, it.col() + size * rspde_order) = it.value();
        }
    }

    // Apply scaling factor and tau
    Q *= std::pow(scale_factor, 2 * alpha) * std::pow(tau, 2);

    return Q;
}

extern "C" {
    void compute_Q_anisotropic(
    double hx, double hy, double hxy, double sigma, double nu,
    const inla_cgeneric_smat_tp *C,
    const inla_cgeneric_smat_tp *Ci,
    const inla_cgeneric_smat_tp *Hxx,
    const inla_cgeneric_smat_tp *Hyy,
    const inla_cgeneric_smat_tp *Hxy,
    double *ret, int rspde_order,
    const inla_cgeneric_mat_tp *rational_table);
}

void compute_Q_anisotropic(double hx, double hy, double hxy, double sigma, double nu,
    const inla_cgeneric_smat_tp *C,
    const inla_cgeneric_smat_tp *Ci,
    const inla_cgeneric_smat_tp *Hxx,
    const inla_cgeneric_smat_tp *Hyy,
    const inla_cgeneric_smat_tp *Hxy,
    const inla_cgeneric_smat_tp *Q_graph,
    double *ret, int rspde_order,
    const inla_cgeneric_mat_tp *rational_table) {
        
        double alpha = nu + 1.0;

        double tau = std::sqrt(std::tgamma(alpha - 1.0) / 
                          (std::tgamma(alpha) * 4.0 * M_PI * std::pow(sigma, 2.0)));
        tau /= std::sqrt(hx * hy * std::sqrt(1.0 - std::pow(hxy, 2.0)));

        // Assembling the FEM matrices
        SparseMatrixColMajor C_eigen = convertInlaToEigen(C);
        SparseMatrixColMajor Ci_eigen = convertInlaToEigen(Ci);
        SparseMatrixColMajor Hxx_eigen = convertInlaToEigen(Hxx);
        SparseMatrixColMajor Hyy_eigen = convertInlaToEigen(Hyy);
        SparseMatrixColMajor Hxy_eigen = convertInlaToEigen(Hxy);
        SparseMatrixColMajor Q_graph_eigen = convertInlaToEigen(Q_graph);
        // Assembling the rational table
        Eigen::MatrixXd rational_table_eigen = inlaToEigenMatrix(rational_table);

        int m_alpha = static_cast<int>(std::floor(alpha));
        
        Q_graph = Q_graph + Eigen::SparseMatrix<double>(Q_graph.transpose());
        
        SparseMatrixColMajor Hxy_transpose = Hxy_eigen.transpose();
        SparseMatrixColMajor L = C_eigen + (hx * hx) * Hxx_eigen + (hy * hy) * Hyy_eigen +
                         (hx * hy * hxy) * (Hxy_eigen + Hxy_transpose);
        SparseMatrixColMajor CiL = Ci_eigen * L;

        SparseMatrixColMajor Q;
        if (std::floor(alpha) == alpha) { // Check if alpha is an integer
            Q = L;
                if (alpha > 1) {
                    for (int k = 1; k < alpha; ++k) {
                        Q = Q * CiL;
                    }
                }
            Q *= tau * tau;
        } else if (rspde_order > 0) {
            Q = anisotropic_precision(L, tau, C_eigen, Ci_eigen, CiL, alpha, m_alpha, rspde_order, 1.0, rational_table_eigen);
        } else {
            throw std::invalid_argument("rspde_order > 0 required");
        }

        Q = Q + 0 * Q_graph_eigen;

        // Extract the values from Q into the result array, using only lower triangular part
        Eigen::SparseMatrix<double, Eigen::ColMajor> Q_triang = Q.triangularView<Eigen::Lower>();

        int count = 0;
        for (int k = 0; k < Q_triang.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(Q_triang, k); it; ++it) {
                ret[count] = it.value();
                count++;
            }
        }        
    }