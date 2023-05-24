#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

extern "C" double logmultnormvdens(int npar, double *entries_mean, 
                        double *entries_prec,
                        double *entries_val);

double logmultnormvdens(int npar, double *entries_mean, 
                        double *entries_prec,
                        double *entries_val) {

                    int i, j, k;
                    
                    Eigen::MatrixXd prec_matrix(npar, npar); 

                        for(i = 0; i < npar; i++){
                            for(j = 0; j < npar; j++){
                                prec_matrix(i,j) = entries_prec[i*npar + j];
                            }
                        }

                    Eigen::VectorXd mean_vec(npar);
                        for(k = 0; k < npar; k++){
                            mean_vec(k) = entries_mean[k];
                        }

                    Eigen::VectorXd val_vec(npar);
                        for(k = 0; k < npar; k++){
                            val_vec(k) = entries_val[k];
                        }

                    Eigen::LLT<Eigen::MatrixXd> chol(prec_matrix);

                    double logdens;

                    Eigen::VectorXd centered_vec(npar);

                    centered_vec = val_vec - mean_vec;

                    logdens = -0.5 * centered_vec.cwiseProduct(prec_matrix * centered_vec).sum();

                    logdens -= npar/2.0 * log(2 * M_PI);

                    logdens += chol.matrixL().toDenseMatrix().diagonal().array().log().sum();

                    return logdens;

                        }