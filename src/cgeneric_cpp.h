#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrixColMajor;

SparseMatrixColMajor convertInlaToEigen(const inla_cgeneric_smat_tp* inlaMat);
Eigen::MatrixXd inlaToEigenMatrix(const inla_cgeneric_mat_tp *mat);
