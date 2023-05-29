#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

extern "C" void compute_Q(int size, double *entries_C, int *i_C, int *j_C,
                        int n_nonzero_C,
                        double *entries_G, int *i_G, int *j_G,
                        int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, int rspde_order, double *theta_entries,
                    double *rat_p, double *rat_r, double rat_k,
                    double *Q_out,
                    int *graph_i, int *graph_j, int M,
                    int matern_par, double start_nu, double nu, double d);

void compute_Q(int size, double *entries_C, int *i_C, int *j_C,
                    int n_nonzero_C,
                    double *entries_G, int *i_G, int *j_G,
                    int n_nonzero_G,
                    double *entries_B_kappa, double *entries_B_tau,
                    int ncol_B, int rspde_order, double *theta_entries,
                    double *rat_p, double *rat_r, double rat_k,
                    double *Q_out,
                    int *graph_i, int *graph_j, int M, 
                    int matern_par, double start_nu, double nu, double d) {

                        double alpha = nu + d/2.0;
                        int m_alpha = (int) floor(alpha);
                        
                        
                        typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_C, trp_G, trp_Q;
                        int k, i, j;

                        
                        // Assemble C and G
                        Eigen::SparseMatrix<double> C(size,size), G(size,size), Q_graph(size*(rspde_order+1), size*(rspde_order+1));

                        for(k = 0; k < n_nonzero_C; k++){
                                trp_C.push_back(Trip(i_C[k],j_C[k],entries_C[k]));
                        }

                        for(k = 0; k < n_nonzero_G; k++){
                                trp_G.push_back(Trip(i_G[k],j_G[k],entries_G[k]));
                        }

                        C.setFromTriplets(trp_C.begin(), trp_C.end());
                        G.setFromTriplets(trp_G.begin(), trp_G.end());                      

                        for(k = 0; k < M; k++){
                                trp_Q.push_back(Trip(graph_i[k],graph_j[k],1));
                        }

                        Q_graph.setFromTriplets(trp_Q.begin(), trp_Q.end());

                        Q_graph = Q_graph + Eigen::SparseMatrix<double>(Q_graph.transpose());

                        // Assemble B_kappa and B_tau

                        Eigen::MatrixXd B_kappa(size, ncol_B), B_tau(size, ncol_B); 

                        for(i = 0; i < size; i++){
                            for(j = 0; j < ncol_B; j++){
                                B_tau(i,j) = entries_B_tau[i*ncol_B + j];
                                B_kappa(i,j) = entries_B_kappa[i*ncol_B + j];
                            }
                        }  

                        if(matern_par == 1){
                            B_kappa.col(0) += 0.5 * log( 8 * nu) * Eigen::VectorXd::Constant(B_kappa.rows(), 1) - 
                            0.5 * log(8 * start_nu) * Eigen::VectorXd::Constant(B_kappa.rows(), 1);
                            B_tau.col(0) += 0.5 * (lgamma(start_nu + d/2.0) -
                                                    lgamma(start_nu) - lgamma(nu + d/2.0) + lgamma(nu)) * Eigen::VectorXd::Constant(B_kappa.rows(),1) +
                                                     start_nu * B_kappa.col(0) - nu * B_kappa.col(0);
                            for(i = 1; i < B_tau.cols(); i++){
                                B_tau.col(i) += start_nu * B_kappa.col(i) - nu * B_kappa.col(i);
                            }
                        } else if(matern_par == 2){
                            B_tau.col(0) += 0.5 * (lgamma(start_nu + d/2.0) -
                                                    lgamma(start_nu) - lgamma(nu + d/2.0) + lgamma(nu)) * Eigen::VectorXd::Constant(B_kappa.rows(),1) +
                                                     start_nu * B_kappa.col(0) - nu * B_kappa.col(0);
                            for(i = 1; i < B_tau.cols(); i++){
                                B_tau.col(i) += start_nu * B_kappa.col(i) - nu * B_kappa.col(i);
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

                        // Scaling L

                        double factor = pow(kappa.minCoeff(),2);

                        L = L / factor;
                        
                        if(m_alpha > 0){
                            CinvL = C.cwiseInverse() * L;
                        }
                        
                        int m;

                        // Assemble first part of Q

                        Eigen::SparseMatrix<double> tau_matrix(size, size);
                        tau_matrix = tau.asDiagonal();

                        Eigen::SparseMatrix<double> Q_tmp(size,size), Q((rspde_order+1)*size, (rspde_order+1)*size);

                        for(k = 0; k < rspde_order; k++){
                            Q_tmp = (L - rat_p[k] * C)/rat_r[k];
                            if(m_alpha>0){
                                for(m = 0; m < m_alpha; m++){
                                    Q_tmp = Q_tmp * CinvL;
                                }
                            }

                            Q_tmp = tau_matrix * Q_tmp * tau_matrix;

                            for (m=0; m < Q_tmp.outerSize(); ++m)
                            {
                                for (Eigen::SparseMatrix<double>::InnerIterator it(Q_tmp,m); it; ++it)
                                {
                                      Q.insert(it.row() + size*k, it.col() + size*k) = it.value();   
                                }
                            }
                        }
                        

                        // Assemble the K part

                        if(m_alpha == 0){
                            Q_tmp = C.cwiseInverse();
                        } else if(m_alpha == 1){
                           Q_tmp = L;      
                        } else{
                            Q_tmp = L;
                            for(m = 0; m < m_alpha-1; m++){
                                Q_tmp = Q_tmp * CinvL;
                            }                               
                        }

                        Q_tmp = tau_matrix * Q_tmp * tau_matrix;


                        for (m=0; m < Q_tmp.outerSize(); ++m) {
                                    for (Eigen::SparseMatrix<double>::InnerIterator it(Q_tmp,m); it; ++it)
                                    {
                                          Q.insert(it.row() + size*rspde_order, it.col() + size*rspde_order) = it.value()/rat_k;   
                                    }
                                }

                       
                        Q = Q + 0 * Q_graph;

                        Q = pow(factor, alpha) * Q;

                        Eigen::SparseMatrix<double> Q_triang((rspde_order+1)*size, (rspde_order+1)*size);
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