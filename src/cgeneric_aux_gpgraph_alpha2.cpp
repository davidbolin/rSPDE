#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

double r_2(double D, double kappa, double tau, int deriv){
    double aD = abs(D);
    double c =  1/(4 * pow(kappa,3) * pow(tau,2));
    double R0 = exp(-kappa * aD);
    if(deriv == 0){
        return(c * (1 + kappa * aD) * R0);
    } else if(deriv == 1){
        return(-pow(kappa,2) * c * D * R0);
    } else{
        return(pow(kappa,2) * c * ( kappa* aD - 1) * R0);
    }
}

extern "C" void compute_Q_alpha2(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int *lower_edges,
                                        int *upper_edges, int lower_edges_len, int upper_edges_len);

void compute_Q_alpha2(int *i_Tc, int *j_Tc, double *x_Tc, double kappa, double tau, int nE, double w,
                                        int nrow_Tc, int ncol_Tc, int n_nonzero_Tc, double *edge_lengths, double *Q_out, int *lower_edges,
                                        int *upper_edges, int lower_edges_len, int upper_edges_len) {
                    
                        typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_Tc;
                        int k, i;

                        
                        // Assemble Tc
                        Eigen::SparseMatrix<double> Tc(nrow_Tc,ncol_Tc);

                        for(k = 0; k < n_nonzero_Tc; k++){
                                trp_Tc.push_back(Trip(i_Tc[k],j_Tc[k],x_Tc[k]));
                        }

                        Tc.setFromTriplets(trp_Tc.begin(), trp_Tc.end());           

                        // Creating auxiliary matrices

                        Eigen::MatrixXd R_00(2,2);

                        int deriv;
                        deriv = 0;
                        R_00(0,0) = r_2(0, kappa = kappa, tau = tau, deriv);
                        deriv = 1;
                        R_00(0,1) = -r_2(0, kappa = kappa, tau = tau, deriv);
                        R_00(1,0) = -r_2(0, kappa = kappa, tau = tau, deriv);
                        deriv = 2;
                        R_00(1,1) = -r_2(0, kappa = kappa, tau = tau, deriv);
                        
                        Eigen::MatrixXd R00i = R_00.inverse();

                        Eigen::MatrixXd R_node = Eigen::MatrixXd::Zero(4,4);
                        
                        R_node.block(0,0,2,2) = R_00;
                        R_node.block(2,2,2,2) = R_00;

                        Eigen::MatrixXd Ajd = Eigen::MatrixXd::Zero(4,4);            

                        Ajd.block(0,0,2,2) = -w * R00i;

                        Ajd.block(2,2,2,2) = -(1-w)*R00i;     

                        // Creating the triplets for Q

                        Eigen::VectorXd i_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd j_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd x_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);

                        int count = 0;
                        double l_e, r_0l, r_11;

                        for(int i=0; i<nE; i++){
                            l_e = edge_lengths[i];
                            deriv = 0;
                            r_0l = r_2(l_e, kappa, tau, deriv);
                            deriv = 2;
                            r_11 = - r_2(l_e, kappa, tau, deriv);

                            Eigen::MatrixXd R_01(2,2);

                            R_01(0,0) = r_0l;
                            deriv = 1;
                            R_01(1,0) = r_2(-l_e, kappa, tau, deriv);
                            R_01(0,1) = r_2(l_e, kappa, tau, deriv);
                            R_01(1,1) = r_11;
                            R_node.block(0,2,2,2) = R_01;
                            R_node.block(2,0,2,2) = R_01.transpose();
                            Eigen::MatrixXd Q_adj = R_node.inverse() + Ajd;


                            // lower edge precision u
                            i_[count] = 4 * i;
                            j_[count] = 4 * i;
                            x_[count] = Q_adj(0, 0);

                            // lower edge  u'
                            i_[count + 1] = 4 * i + 1;
                            j_[count + 1] = 4 * i + 1;
                            x_[count + 1] = Q_adj(1, 1);

                            // upper edge  u
                            i_[count + 2] = 4 * i + 2;
                            j_[count + 2] = 4 * i + 2;
                            x_[count + 2] = Q_adj(2, 2);

                            // upper edge  u'
                            i_[count + 3] = 4 * i + 3;
                            j_[count + 3] = 4 * i + 3;
                            x_[count + 3] = Q_adj(3, 3);

                            // lower edge  u, u'
                            i_[count + 4] = 4 * i;
                            j_[count + 4] = 4 * i + 1;
                            x_[count + 4] = Q_adj(0, 1);
                            i_[count + 5] = 4 * i + 1;
                            j_[count + 5] = 4 * i;
                            x_[count + 5] = Q_adj(0, 1);

                            // upper edge  u, u'
                            i_[count + 6] = 4 * i + 2;
                            j_[count + 6] = 4 * i + 3;
                            x_[count + 6] = Q_adj(2, 3);
                            i_[count + 7] = 4 * i + 3;
                            j_[count + 7] = 4 * i + 2;
                            x_[count + 7] = Q_adj(2, 3);

                            // lower edge  u, upper edge  u,
                            i_[count + 8]  = 4 * i;
                            j_[count + 8]  = 4 * i + 2;
                            x_[count + 8]  = Q_adj(0, 2);
                            i_[count + 9] = 4 * i + 2;
                            j_[count + 9] = 4 * i;
                            x_[count + 9] = Q_adj(0, 2);

                            // lower edge  u, upper edge  u',
                            i_[count + 10] = 4 * i;
                            j_[count + 10] = 4 * i + 3;
                            x_[count + 10] = Q_adj(0, 3);
                            i_[count + 11] = 4 * i + 3;
                            j_[count + 11] = 4 * i;
                            x_[count + 11] = Q_adj(0, 3);

                            // lower edge  u', upper edge  u,
                            i_[count + 12] = 4 * i + 1;
                            j_[count + 12] = 4 * i + 2;
                            x_[count + 12] = Q_adj(1, 2);
                            i_[count + 13] = 4 * i + 2;
                            j_[count + 13] = 4 * i + 1;
                            x_[count + 13] = Q_adj(1, 2);

                            // lower edge  u', upper edge  u',
                            i_[count + 14] = 4 * i + 1;
                            j_[count + 14] = 4 * i + 3;
                            x_[count + 14] = Q_adj(1, 3);
                            i_[count + 15] = 4 * i + 3;
                            j_[count + 15] = 4 * i + 1;
                            x_[count + 15] = Q_adj(1, 3);

                            count += 16;
                        }

                        if(lower_edges_len > 0){
                            for(i = 0; i<lower_edges_len; i++){
                                int ind1 = 4 * (lower_edges[i]-1);
                                int ind2 = 4 * (lower_edges[i]-1) + 1;
                                double x1 = w / R_00(0,0);
                                double x2 = w/R_00(1,1);
                                i_[count] = ind1;
                                i_[count + 1] = ind2;
                                j_[count] = ind1;
                                j_[count+1] = ind2;
                                x_[count] = x1;
                                x_[count+1] = x2;
                                count += 2;
                            }
                        }

                        if(upper_edges_len > 0){
                            for(i = 0; i<upper_edges_len; i++){
                                int ind1 = 4 * (upper_edges[i]-1) + 2;
                                int ind2 = 4 * (upper_edges[i]-1) + 3;
                                double x1 = (1-w) / R_00(0,0);
                                double x2 = (1-w)/R_00(1,1);
                                i_[count] = ind1;
                                i_[count + 1] = ind2;
                                j_[count] = ind1;
                                j_[count+1] = ind2;
                                x_[count] = x1;
                                x_[count+1] = x2;
                                count += 2;
                            }
                        }                 

                        std::vector<Trip> trp_Q;                               

                        Eigen::SparseMatrix<double> Q(4*nE,4*nE);


                        for(k = 0; k < count; k++){
                                trp_Q.push_back(Trip(i_[k],j_[k],x_[k]));
                        }

                        Q.setFromTriplets(trp_Q.begin(), trp_Q.end());        

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