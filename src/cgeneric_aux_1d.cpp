#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <cmath>
#include <numbers>
#include "cgeneric_defs.h"
#include <iostream>

double matern_covariance(double h, double kappa, double nu, double sigma) {
    double out;
    if(h == 0) {
        out = pow(sigma,2.0);
    } else {
        if (nu == 0.5) {
            out = pow(sigma,2.0) * exp(-kappa * abs(h));
        } else {
            out = (pow(sigma,2.0) / (pow(2.0, nu-1.0) * tgamma(nu))) *
                pow(kappa * abs(h), nu);
            double tmp = std::cyl_bessel_k(abs(nu), kappa * abs(h));
            out *= tmp;
        }    
    }
    return(out);
}

double matern_derivative(double h, double kappa, double nu, double sigma, int deriv) 
{
    double out;
    double C = -pow(kappa,2.0) /(2.0 * (nu - 1.0));
    if(deriv == 0) {
        out = matern_covariance(h, kappa, nu, sigma);
    } else if (deriv == 1) {
        if(h==0) {
            out = 0.0;
        } else {
            out = C * h * matern_covariance(h, kappa, nu - 1.0, sigma);
        }
    }
    else if (deriv == 2) {
        out = matern_covariance(h, kappa, nu - 1.0, sigma) + 
            h * matern_derivative(h, kappa, nu - 1.0, sigma, 1);
        out *= C;
    } else {
        out = (deriv - 1.0) * matern_derivative(h, kappa, nu - 1.0, sigma, deriv - 2) + 
                 h * matern_derivative(h, kappa, nu - 1.0, sigma, deriv - 1);
        out *= C;
    }
    return(out);
}


double matern_p(double s,double t,double kappa, double p,double alpha){
    double h = s-t;
    double out, sigma;
    if(p==0){
        out = matern_covariance(h, kappa, alpha - 0.5, 1.0);
    } else {
        double ca = tgamma(alpha)/tgamma(alpha-0.5);
        int fa = floor(alpha);
        double kp = kappa*sqrt(1.0-p);
        out = matern_covariance(h, kp, 0.5, sqrt(ca*sqrt(M_PI)/sqrt(1.0-p)));
        if(alpha >= 1) {
            for(int j =1; j<=fa; j++) {
                sigma = sqrt(ca*tgamma(j-0.5)/tgamma(j));
                out -= matern_covariance(h, kappa, j-0.5, sigma)/pow(p,1.0 - j);
            }
            out = out/pow(p,fa);
        }
    }
    return(out);
}

double matern_p_deriv(double s,double t,double kappa,double p,double alpha,double deriv){
    double h = s-t;
    double out, sigma;
    if(deriv ==0){
        out = matern_p(s,t,kappa,p,alpha);
    } else {
        if(p==0){
            out = matern_derivative(h, kappa, alpha - 0.5, 1, deriv);
        } else {
            double ca = tgamma(alpha)/tgamma(alpha-0.5);
            int fa = floor(alpha);
            double kp = kappa*sqrt(1.0-p);
            out = matern_derivative(h, kp, 0.5, sqrt(ca*sqrt(M_PI)/sqrt(1.0-p)), deriv);
            if(alpha >= 1) {
                out = out/pow(p,fa);
                for(int j = 1; j<= fa; j++) {
                    sigma = sqrt(ca*tgamma(j-0.5)/tgamma(j));
                    out -= matern_derivative(h, kappa, j-0.5, sigma, deriv)/pow(p,fa + 1 - j);
                }
            }
        } 
    }
    return(out);
}

Eigen::MatrixXd matern_p_joint(double s,double t,double kappa, double p, double alpha){
    
    int fa;
    double tmp;
    if(floor(alpha) == alpha) {
        fa = alpha;
    } else {
        fa = floor(alpha) + 1;
    }
    Eigen::MatrixXd mat(fa,fa);
    for(int i = 0; i < fa; i++) {
        for(int j =i; j< fa; j++) {
            if(i==j) {
                mat(i,i) = pow(-1.0,i)*matern_p_deriv(s,t,kappa,p,alpha, 2*i);
            } else {
                tmp = matern_p_deriv(s,t,kappa,p,alpha, i + j);
                mat(i,j) = pow(-1.0,j)*tmp;
                mat(j,i) = pow(-1.0,i)*tmp;
            }
        }
    }
    return(mat);
}

double matern_k(double s,double t,double kappa,double alpha){
    double h = s-t;
    double ca = tgamma(alpha)/tgamma(alpha-0.5);
    double out;
    if(alpha<1){
        if(h==0) {
            out = sqrt(4*M_PI)*ca/kappa;
        } else {
            out = 0;
        }
    } else {
        double fa = floor(alpha);
        double cfa = tgamma(fa)/tgamma(fa-0.5);
        out = matern_covariance(h, kappa, fa - 1/2, sqrt(ca/cfa));
    }
    return(out);
}

double matern_k_deriv(double s,double t,double kappa, double alpha, int deriv){
    double h = s-t;
    double ca = tgamma(alpha)/tgamma(alpha-0.5);
    double out;
    if(alpha<1){
        if(h==0){
            out = sqrt(4*M_PI)*ca/kappa;
        } else {
            out = 0;
        }
    } else {
        int fa = floor(alpha);
        double cfa = tgamma(fa)/tgamma(fa-0.5);
        out = matern_derivative(h, kappa, fa - 0.5, sqrt(ca/cfa), deriv);
    }
    return(out);
}

Eigen::MatrixXd matern_k_joint(double s,double t,double kappa,double alpha){
    
    int fa = floor(alpha);
    int size;
    if(fa == 0) {
        size = 1;
    } else {
        size = fa;
    }
    Eigen::MatrixXd mat(size,size);
    if(fa == 0) {
        mat(0,0) = matern_k(s,t,kappa,alpha);
    } else {
        for(int i = 0; i< fa; i++) {
            for(int j=i; j< fa; j++) {
                if(i==j) {
                    mat(i,i) = pow(-1.0,i)*matern_k_deriv(s,t,kappa,alpha, 2*i);
                } else {
                    mat(i,j) = pow(-1.0,j)*matern_k_deriv(s,t,kappa,alpha, i + j);
                    mat(j,i) = pow(-1.0,i)*matern_k_deriv(s,t,kappa,alpha, i + j);
                }
            }
        }
    }
    return(mat);
}


std::tuple< Eigen::VectorXd, Eigen::VectorXd, 
            Eigen::VectorXd, Eigen::VectorXd> matern_k_chol(int n, double *loc,
                                                            double kappa, 
                                                            int equally_spaced, 
                                                            double alpha) {
    int fa = floor(alpha);
    int N;
    if(fa == 0) {
        fa = 1;
    } 
    N = n*pow(fa,2) + (n-1)*pow(fa,2) - n*fa*(fa -1)/2;
    Eigen::VectorXd ii(N);
    Eigen::VectorXd jj(N);
    Eigen::VectorXd val(N);
    Eigen::VectorXd Fval(fa*n);

    Eigen::MatrixXd Sigma(2*fa,2*fa);
    
    Eigen::MatrixXd Stransp(fa,fa);
    for(int i=0;i<fa;i++){
        for(int j=0;j<fa;j++){
            Stransp(i,j) = pow(-1,i+j);
        }
    }
    //std::cout << "fa = " << fa << ", n = " << n << std::endl;
    Eigen::MatrixXd Sdiag = matern_k_joint(0,0,kappa,alpha);
    //std::cout << "Sdiag = " << std::endl << Sdiag  << std::endl;
    Sigma.topLeftCorner(fa,fa) = Sdiag;
    Sigma.bottomRightCorner(fa,fa) = Sdiag;
    Eigen::MatrixXd Sigma12(fa,fa);
    
    if(equally_spaced==1){
        Sigma12 = matern_k_joint(loc[0],loc[1],kappa,alpha);
        Sigma.topRightCorner(fa,fa) = Sigma12;
        Sigma12 = matern_k_joint(loc[1],loc[0],kappa,alpha);
        Sigma.bottomLeftCorner(fa,fa) = Stransp.cwiseProduct(Sigma12);
    }
    //std::cout << "Sigma = " << std::endl << Sigma  << std::endl;
    int counter = 1;
    int counter2 = 0;
    Eigen::VectorXd tmp;
    for(int i = 0; i < n; i++){
        //std::cout << "i = " << i << std::endl;
        if(i==0){
            Fval(0) = 1/Sdiag(0,0);
            val(0) = 1;
            ii(0) = 0;
            jj(0) = 0;
            if(fa > 1) {
                for(int k = 2; k <= fa; k++) {
                    tmp = Sdiag.topLeftCorner(k-1,k-1).inverse() * Sdiag.block(0,k-1,k-1,1);

                    val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                    val(counter2 + k) = 1;
                    
                    ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                    jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                    double tmpk = Sdiag(k-1,k-1) - Sdiag(Eigen::seq(0,k-2),k-1).dot(tmp);
                    Fval(k) = 1.0/tmpk;
                    counter += 1;
                    counter2 += k;
                }    
            }
        } else {
            if(equally_spaced == 0){
                Sigma12 = matern_k_joint(loc[i-1],loc[i],kappa,alpha);
                Sigma.topRightCorner(fa,fa) = Sigma12;
                Sigma.bottomLeftCorner(fa,fa) = Stransp.cwiseProduct(Sigma12);
            } 
 
            for(int k = fa+1; k <= 2*fa; k++) {
                //std::cout << "k = " << k << ", counter = " << counter << ", counter2 = " << counter2 << std::endl;
                tmp = Sigma.topLeftCorner(k-1,k-1).inverse() * Sigma.block(0,k-1,k-1,1);
                val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                val(counter2 + k) = 1;
                ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                //std::cout << "set tmpk" << std::endl;
                //Fval(counter) = 1.0/(Sigma(k-1,k-1) - Sigma(Eigen::seq(0,k-2),k-1).dot(tmp));
                double tmpk = Sigma(k-1,k-1) - Sigma(Eigen::seq(0,k-2),k-1).dot(tmp);
                Fval(counter) = 1.0/tmpk;
                counter += 1;
                counter2 += k;
            }
        }
    }
    return {ii, jj, val, Fval};
}

std::tuple< Eigen::VectorXd, Eigen::VectorXd, 
            Eigen::VectorXd, Eigen::VectorXd> matern_p_chol(int n, double *loc,
                                                            double kappa,
                                                            double p, 
                                                            int equally_spaced, 
                                                            double alpha) {
    
    int fa, N; 
    if(floor(alpha) == alpha) {
        fa = alpha;
    } else {
        fa = floor(alpha) + 1;
    }
    N = n*pow(fa,2) + (n-1)*pow(fa,2) - n*fa*(fa -1)/2;
    Eigen::VectorXd ii(N);
    Eigen::VectorXd jj(N);
    Eigen::VectorXd val(N);
    Eigen::VectorXd Fval(fa*n);
    
    Eigen::MatrixXd Sigma(2*fa,2*fa);
    
    Eigen::MatrixXd Stransp(fa,fa);
    for(int i=0;i<fa;i++){
        for(int j=0;j<fa;j++){
            Stransp(i,j) = pow(-1.0,i+j);
        }
    }

    Eigen::MatrixXd Sdiag = matern_p_joint(0,0,kappa,p, alpha);
    Eigen::MatrixXd Sod = matern_p_joint(loc[0],loc[1],kappa,p,alpha);
    double di = abs(loc[1]-loc[0]);
        
    Sigma.topLeftCorner(fa,fa) = Sdiag;
    Sigma.bottomRightCorner(fa,fa) = Sdiag;
    Sigma.topRightCorner(fa,fa) = Sod;
    Sigma.bottomLeftCorner(fa,fa) = Stransp.cwiseProduct(Sod);
    
    int counter = 1;
    int counter2 = 0; 
    Eigen::VectorXd tmp;
    
    for(int i = 0; i < n; i++){
        if(i==0){
            Fval(0) = 1/Sdiag(0,0);
            val(0) = 1;
            ii(0) = 0;
            jj(0) = 0;
            
            if(fa > 1) {
                for(int k = 2; k <= fa; k++) {
                    tmp = Sdiag.topLeftCorner(k-1,k-1).inverse() * Sdiag.block(0,k-1,k-1,1);
                    val(Eigen::seq(counter2 +1,counter2 +k-1)) = -tmp;
                    val(counter2 + k) = 1;
                    ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                    jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);

                    //Fval(k) = 1.0/(Sdiag(k-1,k-1) - Sdiag(Eigen::seq(0,k-2),k-1).dot(tmp));
                    double tmpp = Sdiag(k-1,k-1) - Sdiag(Eigen::seq(0,k-2),k-1).dot(tmp);
                    Fval(counter) = 1.0/tmpp;
                    
                    counter += 1;
                    counter2 += k;
                }  
            }
        } else {
            if(equally_spaced == 0){
                if(!(di == abs(loc[i]-loc[i-1]))) {
                    di = abs(loc[i]-loc[i-1]);
                    Sod = matern_p_joint(loc[i-1],loc[i],kappa,p,alpha);
                    Sigma.topRightCorner(fa,fa) = Sod;
                    Sigma.bottomLeftCorner(fa,fa) = Stransp.cwiseProduct(Sod);
                }
            } 
            for(int k = fa+1; k <= 2*fa; k++) {
                tmp = Sigma.topLeftCorner(k-1,k-1).inverse() * Sigma.block(0,k-1,k-1,1);
                val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                val(counter2 + k) = 1;
                ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                double tmpp = Sigma(k-1,k-1) - Sigma(Eigen::seq(0,k-2),k-1).dot(tmp);
                Fval(counter) = 1.0/tmpp;
                    
                counter += 1;
                counter2 += k;
            }
        }
    }
    return {ii, jj, val, Fval};
}

int map_index(int i, int n, int fa1, int fa2, int fa_base1, int fa_base2) {
    int out;
    if(i < n*fa1) { 
        int obs = floor(i/fa1);
        out = i + obs*(fa_base1 - fa1);
    } else {
        int i2 = i - n*fa1;
        int p =floor(i2/(n*fa_base2));
        int obs = floor((i2 - p*n*fa2)/fa2);
        out = i2 + n*(fa_base1 + p*(fa_base2 - fa2)) + obs*(fa_base2 - fa2);
    }
    return(out);
}

extern "C" void compute_Q1d(int n, double *loc, int rspde_order, double kappa,
                           double sigma, double *rat_p, double *rat_r, double rat_k,
                           double *Q_out, int *graph_i, int *graph_j, double nu, int M,
                           int equally_spaced, double nu_upper_bound, int N, double *lconst);

void compute_Q1d(int n, double *loc, int rspde_order, double kappa,
                 double sigma, double *rat_p, double *rat_r, double rat_k,
                    double *Q_out, int *graph_i, int *graph_j, double nu, int M,
                    int equally_spaced, double nu_upper_bound, int N, double *lconst) 
    {
    
    //std::cout << "kappa = " << kappa << ", sigma = " << sigma << ", nu = " << nu << ", nu_ub = " << nu_upper_bound << std::endl;
    //std::cout << "rat_p = " << *rat_p << ", rat_r = " << *rat_r << ", rat_k = " << rat_k << ", rspde_order = " << rspde_order << std::endl;
    double alpha = nu + 0.5;
    double alpha_ub = nu_upper_bound + 0.5;
    
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_L, trp_D, trp_Q;
    
    // setup mapping 
    int fa = floor(alpha);
    if(fa == 0) {
        fa = 1;
    }
    int fa2;
    if(floor(alpha) == alpha) {
        fa2 = alpha;
    } else {
        fa2 = floor(alpha) + 1;
    }
    int fa_base = floor(nu_upper_bound + 0.5);
    int fa_base2;
    if(floor(alpha_ub) == alpha_ub) {
        fa_base2 = alpha_ub;
    } else {
        fa_base2 = floor(alpha_ub) + 1;
    }
    
    Eigen::VectorXd vec_ind(N); 
    int J = fa*n + rspde_order*fa2*n;
    for(int j = 0; j < J; j++) {
        vec_ind(map_index(j, n, fa, fa2, fa_base, fa_base2)) = 1;
    }
        
    // K part
    //std::cout << "K part" << std::endl;
    int m, k, i;
    
    int N1 = n*pow(fa,2) + (n-1)*pow(fa,2) - n*fa*(fa -1)/2;
    
    auto [ii, jj, val, FvalK] = matern_k_chol(n, loc, kappa, equally_spaced, alpha);
    //std::cout << "K chol done" << std::endl;
    for(int k = 0; k < N1; k++){
        trp_L.push_back(Trip(map_index(ii[k], n, fa, fa2, fa_base, fa_base2),
                             map_index(jj[k], n, fa, fa2, fa_base, fa_base2),
                             val[k]));
    }
    int ind;
    //std::cout << "Set D, rat_k = " << rat_k << ", sigma = " << sigma << std::endl;
    for(int k = 0; k < fa*n; k++){
        ind = map_index(k, n, fa, fa2, fa_base, fa_base2);
        trp_D.push_back(Trip(ind, ind, FvalK[k]/(rat_k*pow(sigma,2.0))));
    }
    lconst[0] = FvalK.array().log().sum() - fa*n*log(rat_k) -2*fa*n*log(sigma);
    
    int ind_start = fa*n; 
    
    // P part
    if(rspde_order > 0) {
        //std::cout << "P part" << std::endl;
        int N2 = n*pow(fa2,2) + (n-1)*pow(fa2,2) - n*fa2*(fa2 -1)/2;
        for(i =0; i < rspde_order; i++) {
            auto [ii2, jj2, val2, Fvalp] = matern_p_chol(n, loc, kappa, rat_p[i],equally_spaced, alpha);
            for(int k = 0; k < N2; k++){
                trp_L.push_back(Trip(map_index(ind_start  + ii2[k], n, fa, fa2, fa_base, fa_base2),
                                     map_index(ind_start + jj2[k], n, fa, fa2, fa_base, fa_base2),
                                     val2[k]));
            }  
            for(int k = 0; k < fa2*n; k++){
                ind = map_index(ind_start + k, n, fa, fa2, fa_base, fa_base2);
                trp_D.push_back(Trip(ind, ind, Fvalp[k]/(rat_r[i]*pow(sigma,2.0))));
            }
            
            lconst[0] += Fvalp.array().log().sum() - fa2*n*log(rat_r[i]) -2*fa2*n*log(sigma);
            ind_start += fa2*n;
        }    
    } 
    
    //std::cout << "Set matrices" << std::endl;
    Eigen::SparseMatrix<double> L(N,N), D(N,N), Q_graph(N,N);
    L.setFromTriplets(trp_L.begin(), trp_L.end());
    D.setFromTriplets(trp_D.begin(), trp_D.end());
    //std::cout << "L = " << std::endl << Eigen::MatrixXd(L) << std::endl;
    //std::cout << "D = " << std::endl << Eigen::MatrixXd(D) << std::endl;
    lconst[0] -= n*(fa + rspde_order*fa2)*log(2.0*M_PI);
    lconst[0] /= 2.0;
    
    Eigen::SparseMatrix<double> Q = L.transpose() * D * L;
    
    //std::cout << "Q = " << std::endl << Eigen::MatrixXd(Q) << std::endl;
    for(k = 0; k < N; k++){
        if(vec_ind(k) == 0) {
            Q.coeffRef(k,k) += 1.0;
        }
    }
    
    //std::cout << "Q = " << std::endl << Eigen::MatrixXd(Q) << std::endl;
    
    for(k = 0; k < M; k++){
        trp_Q.push_back(Trip(graph_i[k],graph_j[k],1));
    }
    Q_graph.setFromTriplets(trp_Q.begin(), trp_Q.end());
    Q_graph = Q_graph + Eigen::SparseMatrix<double>(Q_graph.transpose());
    
    Q = Q + 0 * Q_graph;
    
    //std::cout << "Q = " << std::endl << Eigen::MatrixXd(Q) << std::endl;
    
    Eigen::SparseMatrix<double> Q_triang(N,N);
    Q_triang = Q.triangularView<Eigen::Lower>();
                        
    int count = 0;
                        
    for (m=0; m < Q_triang.outerSize(); ++m) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang,m); it; ++it) {                                                                                                                  
            Q_out[count] = it.value();
            count++;
        }
    }
}



