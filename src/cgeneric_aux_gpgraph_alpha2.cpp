#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

double r_2(double D, double kappa, double tau, int deriv){
    double aD = std::abs(D);
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

// sinh(x) - x, computed stably for all x.
// For small |x|, direct evaluation cancels catastrophically (subtracts x from sinh(x) which differ by O(x^3)).
static double sinh_minus_x(double x) {
    if (std::abs(x) < 1.0) {
        double x2 = x * x;
        double term = x * x2 / 6.0;        // n=0: x^3/6
        double sum = term;
        for (int n = 0; n < 8; ++n) {       // n=1..8: x^{2n+3}/(2n+3)!
            term *= x2 / static_cast<double>((2 * n + 4) * (2 * n + 5));
            sum += term;
        }
        return sum;
    }
    return std::sinh(x) - x;
}

// kl*cosh(kl) - sinh(kl), computed stably for all kl (= -(sinh(kl) - kl*cosh(kl))).
static double kl_cosh_minus_sinh(double kl) {
    if (std::abs(kl) < 1.0) {
        double k2 = kl * kl;
        double term = kl * k2 / 3.0;        // n=1: kl^3/3
        double sum = term;
        for (int n = 1; n <= 8; ++n) {      // n=2..9
            term *= k2 / static_cast<double>(n * 2 * (2 * n + 3));
            sum += term;
        }
        return sum;
    }
    return kl * std::cosh(kl) - std::sinh(kl);
}

// Closed-form 4x4 precision matrix Q_adj for the Whittle-Matern alpha=2 SPDE
// restricted to one edge of length l_e. The ordering of the 4 variables is
// (u(0), u'(0), u(l), u'(l)). This matches the R implementation Qalpha2/Q00
// in MetricGraph and is numerically stable for small kappa*l_e (where
// inverting the 4x4 covariance matrix loses all precision).
static void Q00_alpha2(double l_e, double kappa, double tau, Eigen::MatrixXd& Q) {
    Q.setZero();
    double kl = kappa * l_e;
    double c1 = 2.0 * kappa * kl;

    double s_2kl = std::sinh(2.0 * kl);
    double s_kl  = std::sinh(kl);
    double c_kl  = std::cosh(kl);

    // Stable evaluations near kl = 0.
    double s2kl_minus_2kl = sinh_minus_x(2.0 * kl);          // sinh(2kl) - 2kl
    double kl_ckl_minus_skl = kl_cosh_minus_sinh(kl);        // kl*cosh(kl) - sinh(kl)
    // -2kl^2 + cosh(2kl) - 1 = 2*(sinh(kl) - kl)*(sinh(kl) + kl)
    double sinh_minus_kl = sinh_minus_x(kl);
    double denom = 2.0 * sinh_minus_kl * (s_kl + kl);

    Q(0, 0) = Q(2, 2) = c1 * kappa + kappa * kappa * s_2kl;
    Q(0, 1) = Q(1, 0) = c1 * kl;
    Q(2, 3) = Q(3, 2) = -c1 * kl;
    Q(0, 2) = Q(2, 0) = -(2.0 * kappa * kappa * s_kl + c1 * kappa * c_kl);
    Q(0, 3) = Q(3, 0) = c1 * s_kl;
    Q(1, 2) = Q(2, 1) = -c1 * s_kl;
    Q(1, 1) = Q(3, 3) = s2kl_minus_2kl;
    // Q(1,3) = -2*(sinh(kl) - kl*cosh(kl)) = 2*(kl*cosh(kl) - sinh(kl))
    Q(1, 3) = Q(3, 1) = 2.0 * kl_ckl_minus_skl;

    double C = 2.0 * kappa * tau * tau / denom;
    Q *= C;
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

                        // R_00 is the stationary 2x2 covariance of (u(0), u'(0)). It is diagonal
                        // for the Whittle-Matern alpha=2 covariance: R_00 = diag(c, kappa^2 * c)
                        // with c = 1/(4 * kappa^3 * tau^2).
                        Eigen::MatrixXd R_00(2,2);

                        int deriv;
                        deriv = 0;
                        R_00(0,0) = r_2(0, kappa = kappa, tau = tau, deriv);
                        deriv = 1;
                        R_00(0,1) = -r_2(0, kappa = kappa, tau = tau, deriv);
                        R_00(1,0) = -r_2(0, kappa = kappa, tau = tau, deriv);
                        deriv = 2;
                        R_00(1,1) = -r_2(0, kappa = kappa, tau = tau, deriv);

                        // Creating the triplets for Q

                        Eigen::VectorXd i_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd j_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);
                        Eigen::VectorXd x_ = Eigen::VectorXd::Zero(nE*16 + 2*lower_edges_len + 2*upper_edges_len);

                        int count = 0;
                        double l_e;
                        Eigen::MatrixXd Q_adj(4, 4);

                        for(int i=0; i<nE; i++){
                            l_e = edge_lengths[i];

                            // Closed-form per-edge precision (= R_node^{-1} + Ajd analytically).
                            // Numerically stable for small kappa*l_e, where the previous
                            // matrix-inversion approach lost all precision and produced
                            // wrong-sign diagonals on graphs with short edges (e.g. mesh-augmented
                            // graphs with closely spaced observations).
                            Q00_alpha2(l_e, kappa, tau, Q_adj);


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
