#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrixColMajor;

SparseMatrixColMajor convertInlaToEigen(const inla_cgeneric_smat_tp* inlaMat);
Eigen::MatrixXd inlaToEigenMatrix(const inla_cgeneric_mat_tp *mat);
Eigen::SparseMatrix<double> anisotropic_precision(const Eigen::SparseMatrix<double>& L,  
                                                  double tau,                            
                                                  const Eigen::SparseMatrix<double>& C,  
                                                  const Eigen::SparseMatrix<double>& Ci, 
                                                  const Eigen::SparseMatrix<double>& CiL,
                                                  double alpha,                        
                                                  int m_alpha,                         
                                                  int rspde_order,                     
                                                  double scale_factor,                 
                                                  const Eigen::MatrixXd& rational_table);