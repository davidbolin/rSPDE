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
        out = pow(sigma,2);
    } else {
        if (nu == 1 / 2) {
            out = pow(sigma,2) * exp(-kappa * abs(h));
        } else {
            out = (pow(sigma,2) / (pow(2, nu-1) * tgamma(nu))) *
                pow(kappa * abs(h), nu) * std::cyl_bessel_k(nu, kappa * abs(h));
        }    
    }
    return(out);
}

double matern_derivative(double h, double kappa, double nu, double sigma, int deriv) 
{
    double out;
    double C = -pow(kappa,2) /(2.0 * (nu - 1));
    if(deriv == 0) {
        out = matern_covariance(h, kappa, nu, sigma);
    } else if (deriv == 1) {
        if(h==0) {
            out = 0;
        } else {
            out = C * h * matern_covariance(h, kappa, nu - 1, sigma);
        }
    }
    else if (deriv == 2) {
        out = matern_covariance(h, kappa, nu - 1, sigma) + 
            h * matern_derivative(h, kappa, nu - 1, sigma, 1);
        out *= C;
    } else {
        out = (deriv - 1) * matern_derivative(h, kappa, nu - 1, sigma, deriv - 2) + 
                 h * matern_derivative(h, kappa, nu - 1, sigma, deriv - 1);
        out *= C;
    }
    return(out);
}


double matern_p(double s,double t,double kappa, double p,double alpha){
    double h = s-t;
    double out, sigma;
    if(p==0){
        out = matern_covariance(h, kappa, alpha - 1/2, 1.0);
    } else {
        double ca = tgamma(alpha)/tgamma(alpha-0.5);
        int fa = floor(alpha);
        double kp = kappa*sqrt(1-p);
        out = matern_covariance(h, kp, 1/2, sqrt(ca*sqrt(M_PI)/sqrt(1-p)));
        if(alpha >= 1) {
            for(int j =1; j<=fa; j++) {
                sigma = sqrt(ca*tgamma(j-0.5)/tgamma(j));
                out -= matern_covariance(h, kappa, j-1/2, sigma)/pow(p,1 - j);
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
            out = matern_derivative(h, kappa, alpha - 1/2, 1, deriv);
        } else {
            double ca = tgamma(alpha)/tgamma(alpha-0.5);
            int fa = floor(alpha);
            double kp = kappa*sqrt(1-p);
            out = matern_derivative(h, kp, 1/2, sqrt(ca*sqrt(M_PI)/sqrt(1-p)), deriv);
            if(alpha >= 1) {
                out = out/pow(p,fa);
                for(int j = 1; j<= fa; j++) {
                    sigma = sqrt(ca*tgamma(j-0.5)/tgamma(j));
                    out -= matern_derivative(h, kappa, j-1/2, sigma, deriv)/pow(p,fa + 1 - j);
                }
            }
        } 
    }
    return(out);
}

void matern_p_joint(double s,double t,double kappa, double p, double alpha, Eigen::MatrixXd mat){
    
    int fa;
    double tmp;
    if(floor(alpha) == alpha) {
        fa = alpha;
    } else {
        fa = floor(alpha) + 1;
    }
    
    for(int i = 0; i < fa; i++) {
        for(int j =i; j< fa; j++) {
            if(i==j) {
                mat(i,i) = pow(-1,i-1)*matern_p_deriv(s,t,kappa,p,alpha, 2*(i-1));
            } else {
                tmp = matern_p_deriv(s,t,kappa,p,alpha, i-1 + j - 1);
                mat(i,j) = pow(-1,j-1)*tmp;
                mat(j,i) = pow(-1,i-1)*tmp;
            }
        }
    }
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
        out = (h==0)*sqrt(4*M_PI)*ca/kappa;
    } else {
        int fa = floor(alpha);
        double cfa = tgamma(fa)/tgamma(fa-0.5);
        out = matern_derivative(h, kappa, fa - 1/2, sqrt(ca/cfa), deriv);
    }
    return(out);
}

void matern_k_joint(double s,double t,double kappa,double alpha, Eigen::MatrixXd mat){
    int fa = floor(alpha);
    if(fa == 0) {
        mat(1,1) = matern_k(s,t,kappa,alpha);
    } else {
        for(int i = 0; i< fa; i++) {
            for(int j=i; j< fa; j++) {
                if(i==j) {
                    mat(i,i) = pow(-1,i-1)*matern_k_deriv(s,t,kappa,alpha, 2*(i-1));
                } else {
                    mat(i,j) = pow(-1,j-1)*matern_k_deriv(s,t,kappa,alpha, i-1 + j - 1);
                    mat(j,i) = pow(-1,i-1)*matern_k_deriv(s,t,kappa,alpha, i-1 + j - 1);
                }
            }
        }
    }
}


void matern_k_chol(int n, 
                   double *loc,
                   double kappa, 
                   int equally_spaced, 
                   double alpha, 
                   Eigen::VectorXd ii, 
                   Eigen::VectorXd jj, 
                   Eigen::VectorXd val, 
                   Eigen::VectorXd Fval) {
    
    int fa = floor(alpha);

    Eigen::MatrixXd Sigma(2*fa,2*fa);
    
    Eigen::MatrixXd Stransp(fa,fa);
    for(int i=0;i<fa;i++){
        for(int j=0;j<fa;j++){
            Stransp(i,j) = pow(-1,i+j);
        }
    }
    
    Eigen::MatrixXd Sdiag(fa,fa);
    matern_k_joint(0,0,kappa,alpha, Sdiag);
    
    Sigma.topLeftCorner(fa,fa) = Sdiag;
    Sigma.bottomRightCorner(fa,fa) = Sdiag;
    Eigen::MatrixXd Sigma12(fa,fa);
    
    if(equally_spaced==1){
        matern_k_joint(loc[1],loc[2],kappa,alpha, Sigma12);
        Sigma.topRightCorner(fa,fa) = Sigma12;
        matern_k_joint(loc[2],loc[1],kappa,alpha, Sigma12);
        Sigma.bottomLeftCorner(fa,fa) = Sigma12;
    }
    int counter = 1;
    int counter2 = 0;
    Eigen::VectorXd tmp;
    for(int i = 0; i < n; i++){
        if(i==0){
            Fval(0) = 1/Sdiag(0,0);
            val(0) = 1;
            ii(0) = 1;
            jj(0) = 1;
            if(fa > 1) {
                for(int k = 1; k < fa; k++) {
                    tmp = Sdiag.topLeftCorner(k-1,k-1).inverse() * Sdiag.block(0,k-1,k-1,1);

                    val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                    val(counter2 + k) = 1;
                    
                    ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                    jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                    
                    Fval(k) = 1.0/(Sdiag(k,k) - Sdiag(Eigen::seq(0,k-1),k-1).dot(tmp));
                    counter += 1;
                    counter2 += k;
                }    
            }
        } else {
            if(equally_spaced == 0){
                matern_k_joint(loc[i-1],loc[i],kappa,alpha, Sigma12);
                Sigma.topRightCorner(fa,fa) = Sigma12;
                Sigma.bottomLeftCorner(fa,fa) = Stransp*Sigma12;
            } 
            for(int k = fa+1; k <= 2*fa; k++) {
                tmp = Sigma.topLeftCorner(k-1,k-1).inverse() * Sigma.block(0,k-1,k-1,1);
                val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                val(counter2 + k) = 1;
                ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                
                Fval(counter) = 1.0/(Sigma(k,k) - Sigma(Eigen::seq(0,k-1),k-1).dot(tmp));
                
                counter += 1;
                counter2 += k;
            }
        }
    }
}

void matern_p_chol(int n, double *loc,double kappa,double p, int equally_spaced, 
                   double alpha, Eigen::VectorXd ii, Eigen::VectorXd jj, 
                   Eigen::VectorXd val, Eigen::VectorXd Fval) {
    
    int fa; 
    if(floor(alpha) == alpha) {
        fa = alpha;
    } else {
        fa = floor(alpha) + 1;
    }
    
    Eigen::MatrixXd Sigma(2*fa,2*fa);
    
    Eigen::MatrixXd Stransp(fa,fa);
    for(int i=0;i<fa;i++){
        for(int j=0;j<fa;j++){
            Stransp(i,j) = pow(-1,i+j);
        }
    }
    
    Eigen::MatrixXd Sdiag(fa,fa);
    matern_p_joint(0,0,kappa,p, alpha, Sdiag);
    
    Eigen::MatrixXd Sod(fa,fa);
    matern_p_joint(loc[0],loc[1],kappa,p,alpha, Sod);
    double di = abs(loc[1]-loc[0]);
        
    Sigma.topLeftCorner(fa,fa) = Sdiag;
    Sigma.bottomRightCorner(fa,fa) = Sdiag;
    Sigma.topRightCorner(fa,fa) = Sod;
    Sigma.bottomLeftCorner(fa,fa) = Stransp*Sod;
        
    int counter = 1;
    int counter2 = 0; 
    Eigen::VectorXd tmp;
    for(int i = 0; i < n; i++){
        if(i==0){
            Fval(0) = 1/Sdiag(0,0);
            val(0) = 1;
            ii(0) = 1;
            jj(0) = 1;
            
            if(fa > 1) {
                for(int k = 1; k < fa; k++) {
                    tmp = Sdiag.topLeftCorner(k-1,k-1).inverse() * Sdiag.block(0,k-1,k-1,1);
                                  
                    val(Eigen::seq(counter2 +1,counter2 +k-1)) = -tmp;
                    val(counter2 + k) = 1;
                    ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                    jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);

                    Fval(k) = 1.0/(Sdiag(k,k) - Sdiag(Eigen::seq(0,k-1),k-1).dot(tmp));
                    counter += 1;
                    counter2 += k;
                }  
            }
        } else {
            if(equally_spaced == 0){
                if(!(di == abs(loc[i]-loc[i-1]))) {
                    di = abs(loc[i]-loc[i-1]);
                    matern_p_joint(loc[i-1],loc[i],kappa,p,alpha,Sod);
                    Sigma.topRightCorner(fa,fa) = Sod;
                    Sigma.bottomLeftCorner(fa,fa) = Stransp*Sod;
                }
            } 
            for(int k = fa+1; k <= 2*fa; k++) {
                tmp = Sigma.topLeftCorner(k-1,k-1).inverse() * Sigma.block(0,k-1,k-1,1);
                val(Eigen::seq(counter2 + 1,counter2 + k-1)) = -tmp;
                val(counter2 + k) = 1;
                ii(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter,counter);
                jj(Eigen::seq(counter2 + 1,counter2 + k)) = Eigen::VectorXd::LinSpaced(k,counter-k+1,counter);
                    
                Fval(counter) = 1.0/(Sigma(k,k) - Sigma(Eigen::seq(0,k-1),k-1).dot(tmp));
                    
                counter += 1;
                counter2 += k;
            }
        }
    }
}


extern "C" void compute_Q1d(int n, double *loc, int rspde_order, double kappa,
                           double sigma, double *rat_p, double *rat_r, double rat_k,
                           double *Q_out, int *graph_i, int *graph_j, double nu, int M,
                           int equally_spaced);

void compute_Q1d(int n, double *loc, int rspde_order, double kappa,
                 double sigma, double *rat_p, double *rat_r, double rat_k,
                    double *Q_out, int *graph_i, int *graph_j, double nu, int M,
                    int equally_spaced) 
    {
    printf("Start Q1d\n");
    double alpha = nu + 0.5;
    
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_L, trp_D, trp_Q;
    
    // K part
    int N, m, k, i;
    int fa = floor(alpha);
    if(fa == 0) {
        N = n;
    } else {
        N = n*pow(fa,2) + (n-1)*pow(fa,2) - n*fa*(fa -1)/2;
    }
    Eigen::VectorXd ii(N);
    Eigen::VectorXd jj(N);
    Eigen::VectorXd val(N);
    Eigen::VectorXd FvalK(N);
    matern_k_chol(n, loc, kappa, equally_spaced, alpha, ii, jj, val, FvalK);
    for(int k = 0; k < N; k++){
        trp_L.push_back(Trip(ii[k],jj[k],val[k]));
    }
    int ind_start = fa*n; 
    
    // P part
    if(floor(alpha) == alpha) {
        fa = alpha;
    } else {
        fa = floor(alpha) + 1;
    }
    
    N = n*pow(fa,2) + (n-1)*pow(fa,2) - n*fa*(fa -1)/2;
    Eigen::VectorXd ii2(N);
    Eigen::VectorXd jj2(N);
    Eigen::VectorXd val2(N);
    Eigen::VectorXd Fvalp(N);
    int nQ = ind_start + rspde_order*N;
    Eigen::VectorXd Dval(nQ);
    Dval(Eigen::seqN(0,fa*n)) = FvalK/(rat_k*pow(sigma,2));
    
    for(i =0; i < rspde_order; i ++) {
        matern_p_chol(n, loc, kappa, rat_p[i],equally_spaced, alpha, ii2, jj2, val2, Fvalp);
        for(int k = 0; k < N; k++){
            trp_L.push_back(Trip(ind_start +ii2[k],ind_start + jj2[k],val2[k]));
        }    
        Dval(Eigen::seqN(ind_start,fa*n)) = Fvalp/(rat_r[i]*pow(sigma,2));
        ind_start += fa*n;
    }
    
    Eigen::SparseMatrix<double> L(nQ,nQ), D(nQ,nQ), Q_graph(nQ,nQ);
    L.setFromTriplets(trp_L.begin(), trp_L.end());
    
    for(k = 0; k < nQ; k++){
        trp_D.push_back(Trip(k,k,Dval[k]));
    }
    D.setFromTriplets(trp_D.begin(), trp_D.end());
    
    Eigen::SparseMatrix<double> Q = L.transpose() * D * L;
    
    for(k = 0; k < M; k++){
        trp_Q.push_back(Trip(graph_i[k],graph_j[k],1));
    }
    
    Q_graph.setFromTriplets(trp_Q.begin(), trp_Q.end());
    
    Q_graph = Q_graph + Eigen::SparseMatrix<double>(Q_graph.transpose());
    
    Q = Q + 0 * Q_graph;

    Eigen::SparseMatrix<double> Q_triang(nQ,nQ);
    Q_triang = Q.triangularView<Eigen::Lower>();
                        
    int count = 0;
                        
    for (m=0; m < Q_triang.outerSize(); ++m) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q_triang,m); it; ++it) {                                                                                                                  
            Q_out[count] = it.value();
            count++;
        }
    }
}



