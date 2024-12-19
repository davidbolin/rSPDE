#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include "cgeneric_defs.h"
#include "cgeneric_cpp.h"


Eigen::VectorXd intrinsic_mean(Eigen::SparseMatrix<double> Lq) {
    
    Eigen::SparseMatrix<double> S = Lq.selfadjointView<Eigen::Lower>();
    
    int n = Lq.rows();
    
    for (int i = n - 1; i >= 0; --i) {
        Eigen::SparseMatrix<double>::ReverseInnerIterator Si(S, i);
        for (Eigen::SparseMatrix<double>::ReverseInnerIterator ij(Lq, i); ij; --ij) {
            
            Eigen::SparseMatrix<double>::ReverseInnerIterator iL(Lq, i);
            Eigen::SparseMatrix<double>::ReverseInnerIterator iS(S, ij.row());
            Si.valueRef() = 0.0;
            while (iL.row() > i) {
                while (iS && (iL.row() < iS.row())) {
                    --iS;
                }
                if (iS && (iL.row() == iS.row())) {
                    Si.valueRef() -= iL.value() * iS.value();
                    --iS;
                }
                --iL;
            }
            if (i == ij.row()) {
                Si.valueRef() += 1 / iL.value();
                Si.valueRef() /= iL.value();
            } else {
                Si.valueRef() /= iL.value();
                while (iS.row() > i) {
                    --iS;
                }
                iS.valueRef() = Si.value();
            }
            --Si;
        }
    }
    
    return S.diagonal();
}

extern "C" void compute_Q_fintrinsic(
            double tau, double nu,
            const inla_cgeneric_smat_tp *C,
            const inla_cgeneric_smat_tp *Ci,
            const inla_cgeneric_smat_tp *G,
            const inla_cgeneric_smat_tp *Q_graph,    
            double *ret, int rspde_order,
            const inla_cgeneric_mat_tp *rational_table,
            int est_nu, int d, int compute_Q, 
            int compute_mean, int compute_logdet,
            double*ld_out, double *mean_out,
            double scaling,
            const inla_cgeneric_smat_tp *D);


void compute_Q_fintrinsic( double tau, double nu,
                           const inla_cgeneric_smat_tp *C,
                           const inla_cgeneric_smat_tp *Ci,
                           const inla_cgeneric_smat_tp *G,
                           const inla_cgeneric_smat_tp *Q_graph,
                           double *ret, int rspde_order,
                           const inla_cgeneric_mat_tp *rational_table,
                           int est_nu, int d,
                           int compute_Q, int compute_mean, int compute_logdet,
                           double*ld_out, double *mean_out,
                           double scaling,
                           const inla_cgeneric_smat_tp *D) {
    
    double alpha = nu + d/2;
    //std::cout << "tau = " << tau << ", nu = " << nu << std::endl;

    // Assembling the FEM matrices
    Eigen::SparseMatrix<double> C_eigen = convertInlaToEigen(C);
    Eigen::SparseMatrix<double> Ci_eigen = convertInlaToEigen(Ci);
    Eigen::SparseMatrix<double> G_eigen = convertInlaToEigen(G);
    Eigen::SparseMatrix<double> D_eigen = convertInlaToEigen(D);
    
    // Assembling the rational table
    Eigen::MatrixXd rational_table_eigen = inlaToEigenMatrix(rational_table);
    
    int m_alpha = static_cast<int>(std::floor(alpha));
    
    Eigen::SparseMatrix<double> L = G_eigen / scaling;
    Eigen::SparseMatrix<double> CiL = Ci_eigen * L;
    
    Eigen::SparseMatrix<double> Q;
    
    bool int_alpha = std::floor(alpha) == alpha && est_nu == 0; // Check if alpha is an integer
    if (int_alpha) { 
        Q = L;
        if (alpha > 1) {
            for (int k = 1; k < alpha; ++k) {
                Q = Q * CiL;
            }
        }
        Q *= tau * tau;
    } else if (rspde_order > 0) {
        Eigen::SparseMatrix<double> Q_graph_eigen = convertInlaToEigen(Q_graph);           
        Eigen::SparseMatrix<double> Q_graph_transpose = Q_graph_eigen.transpose(); 
        Q_graph_eigen = Q_graph_eigen + Q_graph_transpose;
        Q = anisotropic_precision(L, tau, C_eigen, Ci_eigen, CiL, alpha, m_alpha, 
                                  rspde_order, scaling, rational_table_eigen);
        Q = Q + 0 * Q_graph_eigen;
    } else {
        throw std::invalid_argument("rspde_order > 0 required");
    }
    
    if(compute_Q == 1) {
        Eigen::SparseMatrix<double> Qadj = Q + tau * tau * D_eigen; // add diagonal correction
        // Extract the values from Q into the result array, using only lower triangular part
        Eigen::SparseMatrix<double> Q_triang = Qadj.triangularView<Eigen::Lower>();
        int count = 0;
        for (int k = 0; k < Q_triang.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang, k); it; ++it) {
                ret[count] = it.value();
                count++;
            }
        }      
    }
    
    if(compute_mean || compute_logdet){
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > R;
        Eigen::SparseMatrix<double> Lq; 
        Eigen::VectorXd Qidiag_perm;
        Eigen::VectorXd Qidiag;
        int n = C_eigen.rows();
        if(int_alpha) {
            // Get submatrix Q[-n,-,n]
            Eigen::SparseMatrix<double> Qsub = Q.topLeftCorner(n-1,n-1);
            R.analyzePattern(Qsub);
            R.factorize(Qsub);
            
            Lq = R.matrixL();
            if(compute_logdet){
                // constant = log(|Q|*/(2pi)^((d-1)/2), |Q|* = d|Qsub|
                //log const = log(|Q|*) - (d-1)log(2pi)/2 
                //          = log(2d) + log(|R|*) - (d-1)log(2pi)/2 
                double ldet = Lq.diagonal().array().log().sum(); 
                ld_out[0] = ldet + log(2*n) - (n - 1) * log(2 * M_PI) / 2;
            }
            if(compute_mean){
                Qidiag_perm = intrinsic_mean(Lq);
                Qidiag = R.permutationPinv() * Qidiag_perm;
                for(int i = 0; i < n - 1; i++){
                    mean_out[i] = Qidiag[i];
                }
                mean_out[n-1] = 0;
            }
        } else {
            int start;
            ld_out[0] = (rspde_order + 1)*log(2*n) - (rspde_order + 1)* (n - rspde_order) * log(2 * M_PI) / 2;
            for(int k = 0; k < rspde_order + 1; k ++) {
                start = k*n;
                Eigen::SparseMatrix<double> Qsub = Q.block(start,start,n-1,n-1);
                R.analyzePattern(Qsub);
                R.factorize(Qsub);
                Lq = R.matrixL();
                if(compute_logdet){
                    ld_out[0] += Lq.diagonal().array().log().sum();
                }
                
                if(compute_mean){
                    Qidiag_perm = intrinsic_mean(Lq);
                    Qidiag = R.permutationPinv() * Qidiag_perm;
                    for(int i = 0; i < n - 1; i++){
                        mean_out[start + i] = Qidiag[i];
                    }
                    mean_out[start + n - 1] = 0;
                }
            }    
        }
    }
}


extern "C" void compute_Q_intrinsic(int size, 
                             double *entries_C, int *i_C, int *j_C, int n_nonzero_C,
                             double *entries_G, int *i_G, int *j_G, int n_nonzero_G,
                             double *theta, double *Q_out, int alpha, 
                             int compute_Q, int compute_mean, int compute_logdet,
                             double *ld_out, double *mean_out);

void compute_Q_intrinsic(int size, 
                         double *entries_C, int *i_C, int *j_C, int n_nonzero_C,
                         double *entries_G, int *i_G, int *j_G, int n_nonzero_G,
                         double *theta, double *Q_out, int alpha, 
                         int compute_Q, int compute_mean, int compute_logdet,
                         double *ld_out, double *mean_out) {
    
    
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_C, trp_G;
    int k, i;
    
    
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
    
        
    double ltau = theta[0];
    double lkappa = theta[1];
        
    double tau2 = exp(2*ltau);
    double kappa2 = exp(2*lkappa);

    // Create vector of the parts of Q
    
    Eigen::SparseMatrix<double> L(size,size), CinvL(size,size);
    
    L = kappa2*C + G;
    
    CinvL = C.cwiseInverse() * L;
    
    Eigen::SparseMatrix<double> Q(size,size);
    
    Q = G*CinvL;
    
    if(alpha ==2){
        Q = Q * CinvL;
    }
    
    Q = tau2 * Q;
    
    if(compute_Q== 1) {
        Eigen::SparseMatrix<double> Q_triang(size, size);
        Q_triang = Q.triangularView<Eigen::Lower>();
        
        
        int count = 0;
        int m;
        for (m=0; m < Q_triang.outerSize(); ++m)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang,m); it; ++it)
            {                                                                                                                  
                Q_out[count] = it.value();
                count++;
            }
        }    
    }
    
    if(compute_mean || compute_logdet){
        // Get submatrix Q[-n,-,n]
        Eigen::SparseMatrix<double> Qsub(size-1,size-1);
        Qsub = Q.topLeftCorner(size-1,size-1);
        
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > R;
        R.analyzePattern(Qsub);
        R.factorize(Qsub);
        Eigen::SparseMatrix<double> L = R.matrixL();
        if(compute_logdet){
            // constant = log(|Q|*/(2pi)^((d-1)/2), |Q|* = d|Qsub|
            //log const = log(|Q|*) - (d-1)log(2pi)/2 
            //          = log(2d) + log(|R|*) - (d-1)log(2pi)/2 
            double ldet = L.diagonal().array().log().sum(); 
            ld_out[0] = ldet + log(2*size) - (size- 1) * log(2 * M_PI) / 2;
        }
        
        
        if(compute_mean){
            Eigen::SparseMatrix<double> S = L.selfadjointView<Eigen::Lower>();
            
            int n = L.rows();
            
            for (int i = n - 1; i >= 0; --i) {
                Eigen::SparseMatrix<double>::ReverseInnerIterator Si(S, i);
                for (Eigen::SparseMatrix<double>::ReverseInnerIterator ij(L, i); ij; --ij) {
                
                    Eigen::SparseMatrix<double>::ReverseInnerIterator iL(L, i);
                    Eigen::SparseMatrix<double>::ReverseInnerIterator iS(S, ij.row());
                    Si.valueRef() = 0.0;
                    while (iL.row() > i) {
                        while (iS && (iL.row() < iS.row())) {
                            --iS;
                        }
                        if (iS && (iL.row() == iS.row())) {
                            Si.valueRef() -= iL.value() * iS.value();
                            --iS;
                        }
                        --iL;
                    }
                    if (i == ij.row()) {
                        Si.valueRef() += 1 / iL.value();
                        Si.valueRef() /= iL.value();
                    } else {
                        Si.valueRef() /= iL.value();
                        while (iS.row() > i) {
                            --iS;
                        }
                        iS.valueRef() = Si.value();
                    }
                    --Si;
                }
            }
            
            Eigen::VectorXd Qidiag_perm = S.diagonal();
            Eigen::VectorXd Qidiag = R.permutationPinv() * Qidiag_perm;
            for(i = 0; i < n; i++){
                mean_out[i] = Qidiag[i];
            }
        }
    }
}




