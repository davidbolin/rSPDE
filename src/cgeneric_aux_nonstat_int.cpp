#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

extern "C" void compute_Q_integer(int size, double *entries_C, int *i_C, int *j_C,
                        int n_nonzero_C,
                        double *entries_G, int *i_G, int *j_G,
                        int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, double *theta_entries,
                    double *Q_out, int alpha);

void compute_Q_integer(int size, double *entries_C, int *i_C, int *j_C,
                        int n_nonzero_C,
                        double *entries_G, int *i_G, int *j_G,
                        int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, double *theta_entries,
                    double *Q_out, int alpha) {
                    
                        
                        typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_C, trp_G;
                        int k, i, j;

                        
                        // Assemble C and G
                        Eigen::SparseMatrix<double> C(size,size), G(size,size);

                        for(k = 0; k < n_nonzero_C; k++){
                                trp_C.push_back(Trip(i_C[k],j_C[k],entries_C[k]));
                        }

                        for(k = 0; k < n_nonzero_G; k++){
                                trp_G.push_back(Trip(i_G[k],j_G[k],entries_G[k]));
                        }

                        C.setFromTriplets(trp_C.begin(), trp_C.end());
                        G.setFromTriplets(trp_G.begin(), trp_G.end());                      

                        // Assemble B_kappa and B_tau

                        Eigen::MatrixXd B_kappa(size, ncol_B), B_tau(size, ncol_B); 

                        for(i = 0; i < size; i++){
                            for(j = 0; j < ncol_B; j++){
                                B_tau(i,j) = entries_B_tau[i*ncol_B + j];
                                B_kappa(i,j) = entries_B_kappa[i*ncol_B + j];
                            }
                        }

                        // get kappa and tau

                        Eigen::VectorXd theta(ncol_B);
                        theta(0) = 1;
                        for(k = 1; k < ncol_B; k++){
                            theta(k) = theta_entries[k-1];
                        }


                        Eigen::VectorXd kappa = (B_kappa * theta).array().exp();
                        Eigen::VectorXd tau = (B_tau * theta).array().exp();
                        
                        // Create vector of the parts of Q

                        Eigen::VectorXd Cdiag = C.diagonal();

                        Eigen::SparseMatrix<double> L(size,size), CinvL(size,size);

                        L = kappa.cwiseProduct(kappa).cwiseProduct(Cdiag).asDiagonal();
                        L = L + G;

                        if(alpha > 1){
                            CinvL = C.cwiseInverse() * L;
                        }
                        
                        int m;

                        // Assemble first part of Q

                        Eigen::SparseMatrix<double> tau_matrix(size, size);
                        tau_matrix = tau.asDiagonal();

                        Eigen::SparseMatrix<double> Q(size,size);

                        Q = L;

                        if(alpha > 1){
                            for(k = 1; k < alpha; k++){
                                    Q = Q * CinvL;
                                }
                        }

                        Q = tau_matrix * Q * tau_matrix;


                        

                        Eigen::SparseMatrix<double> Q_triang(size, size);
                        Q_triang = Q.triangularView<Eigen::Lower>();
                        

                        int count = 0;
                        
                        for (m=0; m < Q_triang.outerSize(); ++m)
                                {
                                    for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang,m); it; ++it)
                                    {                                                                                                                  
                                          Q_out[count] = it.value();
                                          count++;
                                    }
                                }

                        }