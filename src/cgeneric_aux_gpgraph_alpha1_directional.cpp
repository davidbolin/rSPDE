#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

extern "C" void compute_Q_alpha1_directional(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int stat_ind_len, int *stat_ind, int BC);

void compute_Q_alpha1_directional(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int stat_ind_len, int *stat_ind, int BC) {
                    
                        typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_Tc;
                        int k, i;

                        
                        // Assemble Tc
                        Eigen::SparseMatrix<double> Tc(nrow_Tc,ncol_Tc);

                        for(k = 0; k < n_nonzero_Tc; k++){
                                trp_Tc.push_back(Trip(i_Tc[k],j_Tc[k],x_Tc[k]));
                        }

                        Tc.setFromTriplets(trp_Tc.begin(), trp_Tc.end());           

                        // Creating the triplets for Q

                        // if(BC > 0){
                            Eigen::VectorXd i_ = Eigen::VectorXd::Zero(nE*4 + stat_ind_len);
                            Eigen::VectorXd j_ = Eigen::VectorXd::Zero(nE*4 + stat_ind_len);
                            Eigen::VectorXd x_ = Eigen::VectorXd::Zero(nE*4 + stat_ind_len);
                        // } else{
                        //     Eigen::VectorXd i_ = Eigen::VectorXd::Zero(nE*4);
                        //     Eigen::VectorXd j_ = Eigen::VectorXd::Zero(nE*4);
                        //     Eigen::VectorXd x_ = Eigen::VectorXd::Zero(nE*4);
                        // }


                        int count = 0;
                        double l_e, c1, c2, one_m_c2, c_1_upper, c_1_lower, c_2;
                        
                        double factor =  2 * kappa * pow(tau,2);

                        for(int i=0; i<nE; i++){
                            l_e = edge_lengths[i];
                            c1 = exp(-kappa*l_e);
                            c2 = pow(c1,2);
                            one_m_c2 = 1-c2;
                            c_1_upper = w + c2/one_m_c2;
                            c_1_lower = (1-w) + c2/one_m_c2;
                            c_2 = -c1/one_m_c2;

                            i_[count] = 2 * i;
                            j_[count] = 2 * i;
                            x_[count] = c_1_upper;

                            i_[count + 1] = 2 * i + 1;
                            j_[count + 1] = 2 * i + 1;
                            x_[count + 1] = c_1_lower;

                            i_[count + 2] = 2 * i;
                            j_[count + 2] = 2 * i + 1;
                            x_[count + 2] = c_2;

                            i_[count + 3] = 2 * i + 1;
                            j_[count + 3] = 2 * i;
                            x_[count + 3] = c_2;

                            count += 4;
                        }

                        if(BC > 0){
                            for (int ii = 0; ii < stat_ind_len; ii++) {
                                int ind = stat_ind[ii];
                                i_[count] = ind;
                                j_[count] = ind;
                                x_[count] = (1-w); 
                                count = count + 1;
                            }
                        }   

                        std::vector<Trip> trp_Q;                               

                        Eigen::SparseMatrix<double> Q(2*nE,2*nE);


                        for(k = 0; k < count; k++){
                                trp_Q.push_back(Trip(i_[k],j_[k],x_[k]));
                        }

                        Q.setFromTriplets(trp_Q.begin(), trp_Q.end());     
                        Q = Q * factor;   

                        Eigen::SparseMatrix<double> Q_tilde = Tc * Q * Tc.transpose();
                        
                        // Prepare Q to be sent to inla
                        Eigen::SparseMatrix<double> Q_triang(Q_tilde.rows(), Q_tilde.cols());
                        Q_triang = Q_tilde.triangularView<Eigen::Lower>();
                        

                        count = 0;
                        
                        for (int m=0; m < Q_triang.outerSize(); ++m)
                                {
                                    for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang,m); it; ++it)
                                    {                                                                                                                  
                                          Q_out[count] = it.value();
                                          count++;
                                    }
                                }

                        }