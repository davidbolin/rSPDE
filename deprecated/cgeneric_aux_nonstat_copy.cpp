#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

// extern "C" void compute_Q(int size, double *entries_C, int *i_C, int *j_C, 
//                         double *entries_G, int *i_G, int *j_G,
//                         double *entries_B_kappa, double *entries_B_tau,
//                         int ncol_B, int rspde_order, 
//                         double *theta_entries, double *rat_p, double *rat_r, 
//                         double rat_k, int m_alpha,
//                         double *Q_out, int *i_Q, int *j_Q,
//                         int *graph_i, int *graph_j,int M);

void compute_Q(int size, double *entries_C, int *i_C, int *j_C, 
                        double *entries_G, int *i_G, int *j_G,
                        double *entries_B_kappa, double *entries_B_tau,
                        int ncol_B, int rspde_order,
                        double *theta_entries, double *rat_p, double *rat_r, 
                        double rat_k, int m_alpha, int *graph_i, int *graph_j, int M) {

                        typedef Eigen::Triplet<double> Trip;
                        std::vector<Trip> trp_C, trp_G, trp_Q;
                        int k, i, j;

                        
                        // Assemble C and G
                        Eigen::SparseMatrix<double> C(size,size), G(size,size), Q_graph(size*(rspde_order+1), size*(rspde_order+1));

                        for(k = 0; k < size; k++){
                                trp_C.push_back(Trip(i_C[k],j_C[k],entries_C[k]));
                                trp_G.push_back(Trip(i_G[k],j_G[k],entries_G[k]));
                        }

                        C.setFromTriplets(trp_C.begin(), trp_C.end());
                        G.setFromTriplets(trp_G.begin(), trp_G.end());

                        std::cout << C << std::endl;
                        std::cout << G << std::endl;

                        for(k = 0; k < M; k++){
                                trp_Q.push_back(Trip(graph_i[k],graph_j[k],1));
                        }

                        Q_graph.setFromTriplets(trp_Q.begin(), trp_Q.end());

                        Q_graph = Q_graph + Eigen::SparseMatrix<double>(Q_graph.transpose());

                        std::cout << Q_graph << std::endl;           

                        // Assemble B_kappa and B_tau

                        Eigen::MatrixXd B_kappa(size, ncol_B), B_tau(size, ncol_B); 

                        for(i = 0; i < size; i++){
                            for(j = 0; j < ncol_B; j++){
                                B_tau(i,j) = entries_B_tau[i*ncol_B + j];
                                B_kappa(i,j) = entries_B_kappa[i*ncol_B + j];
                            }
                        }

                        std::cout << B_tau << std::endl;

                        std::cout << B_kappa << std::endl;

                        // get kappa and tau

                        Eigen::VectorXd theta(size);
                        for(k = 0; k < size; k++){
                            theta(k) = theta_entries[k];
                        }

                        std::cout << theta << std::endl;

                        Eigen::VectorXd kappa = (B_kappa * theta).array().exp();
                        Eigen::VectorXd tau = (B_tau * theta).array().exp();

                        std::cout << kappa << std::endl;
                        std::cout << tau << std::endl;

                        // Create vector of the parts of Q

                        Eigen::VectorXd Cdiag = C.diagonal();

                        std::cout << Cdiag << std::endl;

                        Eigen::SparseMatrix<double> L(size,size), CinvL(size,size);


                        L = kappa.cwiseProduct(kappa).cwiseProduct(Cdiag).asDiagonal();
                        L = L + G;
                        
                        std::cout << L << std::endl;

                        if(m_alpha > 0){
                            CinvL = C.cwiseInverse() * L;
                        }

                        std::cout << CinvL << std::endl;
                        
                        int m;

                        // Assemble first part of Q

                        Eigen::SparseMatrix<double> Q_tmp(size,size), Q((rspde_order+1)*size, (rspde_order+1)*size);

                        for(k = 0; k < rspde_order; k++){
                            Q_tmp = (L - rat_p[k] * C)/rat_r[k];
                            if(m_alpha>0){
                                for(m = 0; m < m_alpha; m++){
                                    Q_tmp = Q_tmp * CinvL;
                                }
                            }

                            for (int m=0; m < Q_tmp.outerSize(); ++m)
                            {
                                for (Eigen::SparseMatrix<double>::InnerIterator it(Q_tmp,m); it; ++it)
                                {
                                      Q.insert(it.row() + size*k, it.col() + size*k) = it.value();   
                                }
                            }
                        }

                        // Assemble the K part

                        if(m_alpha == 0){
                            Q_tmp = C;
                        } else if(m_alpha == 1){
                           Q_tmp = L;      
                        } else{
                            Q_tmp = L;
                            for(m = 0; m < m_alpha-1; m++){
                                Q_tmp = Q_tmp * CinvL;
                            }                               
                        }


                        for (int m=0; m < Q_tmp.outerSize(); ++m)
                                {
                                    for (Eigen::SparseMatrix<double>::InnerIterator it(Q_tmp,m); it; ++it)
                                    {
                                          Q.insert(it.row() + size*rspde_order, it.col() + size*rspde_order) = it.value()/rat_k;   
                                    }
                                }

                        std::cout << Q << std::endl;
                        
                        // Q = Q + 0 * Q_graph;


                        // for (int m=0; m < Q.outerSize(); ++m)
                        //         {
                        //             for (Eigen::SparseMatrix<double>::InnerIterator it(Q,m); it; ++it)
                        //             {
                        //                   Q_out[it.row()*m + it.col()] = it.value();
                        //             }
                        //         }

                        }

int main(){
    int size = 3;
    double entries_C[3];
    entries_C[0] = 0.5;
    entries_C[1] = 0.5;
    entries_C[2] = 0.5;

    int i_C[3], j_C[3];
    i_C[0] = j_C[0] = 0;
    i_C[1] = j_C[1] = 1;
    i_C[2] = j_C[2] = 2;

    double entries_G[3];
    int i_G[3], j_G[3];
    entries_G[0] = 1;
    entries_G[1] = 1;
    entries_G[2] = 1;
    i_G[0] = j_G[0] = 0;
    i_G[1] = j_G[1] = 1;
    i_G[2] = j_G[2] = 2;

    double entries_B_kappa[9] = {1, 0, 1, 1, 0, 1, 1, 0, 1};
    double entries_B_tau[9] = {1, 1, 0, 1, 1, 0, 1, 1, 0};
    int ncol_B = 3;

    int rspde_order = 2;
    
    double theta_entries[3] = {0, 0, 0};

    double rat_p[3] = {0.5, 0.5, 0.5};

    double rat_r[3] = {1, 1, 1};

    double rat_k = 2;

    int m_alpha = 2;

    int graph_i[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    int graph_j[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    int M = 9;

    compute_Q(size, entries_C, i_C, j_C, 
                        entries_G, i_G, j_G,
                        entries_B_kappa, entries_B_tau,
                        ncol_B, rspde_order,
                        theta_entries, rat_p, rat_r, 
                        rat_k, m_alpha, graph_i, graph_j, M);



    return 0;

}