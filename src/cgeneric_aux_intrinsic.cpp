#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

extern "C" void compute_Q_intrinsic(int size, 
                                 double *entries_C, int *i_C, int *j_C, int n_nonzero_C,
                                 double *entries_G, int *i_G, int *j_G, int n_nonzero_G,
                                 double *theta, double *Q_out, int alpha, 
                                 int compute_Q, int compute_mean, int compute_logdet,
                                 double*ld, double *mean);

void compute_Q_intrinsic(int size, 
                       double *entries_C, int *i_C, int *j_C, int n_nonzero_C,
                       double *entries_G, int *i_G, int *j_G, int n_nonzero_G,
                       double *theta, double *Q_out, int alpha, 
                       int compute_Q, int compute_mean, int compute_logdet,
                       double*ld_out, double *mean_out) {
    
    
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
        
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double, 0, int> > R;
        R.analyzePattern(Qsub);
        R.factorize(Qsub);
        Eigen::SparseMatrix<double, 0, int> L = R.matrixL();
        if(compute_logdet){
            // constant = log(|Q|*/(2pi)^((d-1)/2), |Q|* = d|Qsub|
            //log const = log(|Q|*) - (d-1)log(2pi)/2 
            //          = log(2d) + log(|R|*) - (d-1)log(2pi)/2 
            double ldet = L.diagonal().array().log().sum(); 
            ld_out[0] = ldet + log(2*size) - (size- 1) * log(2 * M_PI) / 2;
        }
        
        
        if(compute_mean){
            Eigen::SparseMatrix<double, 0, int> S = L.selfadjointView<Eigen::Lower>();
            
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




