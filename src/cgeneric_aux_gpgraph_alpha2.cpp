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

double Q00(int i, int j, double l, double kappa, double tau){
    
    double kl = kappa * l;
    double c1 =  2.0 * kappa * kl;
    double C = 2.0 * kappa * pow(tau,2) / (-2.0*pow(kl,2) + cosh(2.0*kl) - 1.0);
    double out = 0;
    if( (i==0 && j == 0) || (i == 2 && j == 2) ){
        out = C * (c1*kappa + pow(kappa,2) * sinh(2.0*kl));
    } else if ( (i == 0 && j == 1) || (i == 1 && j == 0) ) {
        out = C * c1 * kl;
    } else if ( (i == 2 && j == 3) || (i == 3 && j == 2) ) {
        out = -C * c1 * kl;
    } else if ( (i == 0 && j == 2) || (i == 2 && j == 0) ) {
        out =  C * (-(2.0*pow(kappa,2) * sinh(kl) + c1*kappa*cosh(kl)));
    }  else if ( (i == 0 && j == 3) || (i == 3 && j == 0) ) {
        out = C* c1 * sinh(kl);
    } else if ( (i == 1 && j == 2) || (i == 2 && j == 1) ) {
        out = -C * c1 * sinh(kl);
    } else if ( (i == 1 && j == 1) || (i == 3 && j == 3) ) {
        out = C * (sinh(2.0*kl) - 2.0*kl);
    } else if ( (i == 1 && j == 3) || (i == 3 && j == 1) ) {
        out = C * (-2.0*(sinh(kl) - kl*cosh(kl)));
    } 
    return(out);
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
                        
                        
                        int deriv;
                        deriv = 0;
                        double R_00 = r_2(0, kappa = kappa, tau = tau, deriv);
                        deriv = 2;
                        double R_11 = -r_2(0, kappa = kappa, tau = tau, deriv);

                        // Creating the triplets for Q

                        Eigen::VectorXd i_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd j_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd x_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);

                        int count = 0;
                        double l_e, Qij;

                        for(int i=0; i<nE; i++){
                            l_e = edge_lengths[i];
                            
                            // lower edge precision u
                            Qij = Q00(0, 0, l_e, kappa, tau);    
                            i_[count] = 4 * i;
                            j_[count] = 4 * i;
                            x_[count] = Qij;
                        
                            // lower edge  u'
                            Qij = Q00(1, 1, l_e, kappa, tau);    
                            i_[count + 1] = 4 * i + 1;
                            j_[count + 1] = 4 * i + 1;
                            x_[count+1] = Qij;
                        

                            // upper edge  u
                            Qij = Q00(2, 2, l_e, kappa, tau);    
                            i_[count + 2] = 4 * i + 2;
                            j_[count + 2] = 4 * i + 2;
                            x_[count+2] = Qij;

                            // upper edge  u'
                            Qij = Q00(3, 3, l_e, kappa, tau);    
                    
                            i_[count + 3] = 4 * i + 3;
                            j_[count + 3] = 4 * i + 3;
                            x_[count+3] = Qij;
                           
                            // lower edge  u, u'
                            Qij = Q00(0, 1, l_e, kappa, tau);    
                            i_[count + 4] = 4 * i;
                            j_[count + 4] = 4 * i + 1;
                            x_[count + 4] = Qij;
                            i_[count + 5] = 4 * i + 1;
                            j_[count + 5] = 4 * i;
                            x_[count + 5] = Qij;

                            // upper edge  u, u'
                            Qij = Q00(2, 3, l_e, kappa, tau);    
                            
                            i_[count + 6] = 4 * i + 2;
                            j_[count + 6] = 4 * i + 3;
                            x_[count + 6] = Qij;
                            i_[count + 7] = 4 * i + 3;
                            j_[count + 7] = 4 * i + 2;
                            x_[count + 7] = Qij;

                            // lower edge  u, upper edge  u,
                            Qij = Q00(0, 2, l_e, kappa, tau);    
                            i_[count + 8]  = 4 * i;
                            j_[count + 8]  = 4 * i + 2;
                            x_[count + 8]  = Qij;
                            i_[count + 9] = 4 * i + 2;
                            j_[count + 9] = 4 * i;
                            x_[count + 9] = Qij;

                            // lower edge  u, upper edge  u',
                            Qij = Q00(0, 3, l_e, kappa, tau);    
                            i_[count + 10] = 4 * i;
                            j_[count + 10] = 4 * i + 3;
                            x_[count + 10] = Qij;
                            i_[count + 11] = 4 * i + 3;
                            j_[count + 11] = 4 * i;
                            x_[count + 11] = Qij;

                            // lower edge  u', upper edge  u,
                            Qij = Q00(1, 2, l_e, kappa, tau);    
                            i_[count + 12] = 4 * i + 1;
                            j_[count + 12] = 4 * i + 2;
                            x_[count + 12] = Qij;
                            i_[count + 13] = 4 * i + 2;
                            j_[count + 13] = 4 * i + 1;
                            x_[count + 13] = Qij;

                            // lower edge  u', upper edge  u',
                            Qij = Q00(1, 3, l_e, kappa, tau);    
                            
                            i_[count + 14] = 4 * i + 1;
                            j_[count + 14] = 4 * i + 3;
                            x_[count + 14] = Qij;
                            i_[count + 15] = 4 * i + 3;
                            j_[count + 15] = 4 * i + 1;
                            x_[count + 15] = Qij;

                            count += 16;
                        }

                        if(lower_edges_len > 0){
                            for(i = 0; i<lower_edges_len; i++){
                                int ind1 = 4 * (lower_edges[i]-1);
                                int ind2 = 4 * (lower_edges[i]-1) + 1;
                                double x1 = w / R_00;
                                double x2 = w/R_11;
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
                                double x1 = (1-w) / R_00;
                                double x2 = (1-w)/R_11;
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