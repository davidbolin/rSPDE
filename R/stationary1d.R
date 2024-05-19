#' Rational approximation of the Matern covariance
#' 
#' Computes a rational approximation of the Matern covariance function on intervals. 
#' @param h Distances to compute the covariance for
#' @param order The order of the approximation
#' @param kappa Range parameter
#' @param nu Smoothness parameter
#' @param sigma Standard deviation
#' @param type_rational Method used to compute the coefficients of the rational approximation.
#' @param type_interp Interpolation method for the rational coefficients. 
#'
#' @return The covariance matrix of the approximation
#' @export
#' @examples
#' h <- seq(from = 0, to = 1, length.out = 100)
#' cov.true <- matern.covariance(h, kappa = 10, sigma = 1, nu = 0.8)
#' cov.approx <- matern.rational(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)
#' 
#' plot(h, cov.true)
#' lines(h, cov.approx, col = 2)
#'  
matern.rational = function(h, 
                           order,
                           kappa,
                           nu, 
                           sigma,
                           type_rational = "brasil", 
                           type_interp = "linear")
{
    
    if(is.matrix(loc) && min(dim(h)) > 1) {
        stop("Only one dimensional locations supported.")
    }
    alpha = nu+1/2
    coeff <- interp_rational_coefficients(order = order, 
                                          type_rational_approx = type_rational, 
                                          type_interp = type_interp, 
                                          alpha = alpha)
    r <- coeff$r
    p <- coeff$p
    k <- coeff$k
    if (nu>0 && nu<0.5)
    {
        sigma_rational = 0
        for (i in 1:length(p)){
            sigma_rational = sigma_rational+ r[i]*sigma^2*kappa^2*exp(-sqrt(kappa^2*(1-p[i]))*abs(h))/(2*sqrt(kappa^2*(1-p[i])))
        }
        sigma_rational[h==0] = sigma_rational[h==0] + k*sigma^2
        sigma_rational = sigma_rational*kappa^(-2*alpha)
    }
    
    else {
        a =floor(alpha)
        b = gamma(a)*sqrt(4*pi)*kappa^(2*a-1)*2^(a-(3/2))
        
        rho1 = ((k*sigma^2*kappa^(2*a))/b)*(kappa*abs(h))^(a-1/2)*besselK(kappa * abs(h), a-1/2)
        rho1[h == 0] <- k*sigma^2*kappa^(2*a)*gamma(a-1/2)/(gamma(a)*sqrt(4*pi)*kappa^(2*a-1))
        tmp = 0
        for (i in 1:length(p))
        {
            rho2= r[i]*sigma^2*kappa^(2*a+2)*exp(-sqrt(kappa^2*(1-p[i]))*abs(h))/((2*((kappa^2*p[i])^(a)))*sqrt(kappa^2*(1-p[i])))
            
            sum=0
            for (j in 1: floor(alpha))
            {
                rho3= (r[i]*sigma^2*kappa^(2*a+2)/((kappa^2*p[i])^(a+1-j)*gamma(j)*sqrt(4*pi)*kappa^(2*j-1)*2^(j-3/2)))*(kappa*abs(h))^(j-1/2)*besselK(kappa * abs(h), j-1/2)
                rho3[h == 0] <- r[i]*sigma^2*kappa^(2*a+2)*gamma(j-1/2)/((kappa^2*p[i])^(a+1-j)*gamma(j)*sqrt(4*pi)*kappa^(2*j-1))
                sum=sum+rho3
            }
            tmp=tmp+(rho2-sum)
        }
        sigma_rational=(rho1+tmp)*kappa^(-2*alpha)
    }
    scaling <- ((gamma(alpha)*sqrt(4*pi)*kappa^(2*nu))/gamma(nu))
    return(scaling*as.matrix(sigma_rational))
}


#' Precision matrix of rational approximation of Matern covariance
#'
#' Computes the precision matrix for a rational approximation of the Matern covariance on intervals.
#' @param loc Locations at which the precision is evaluated
#' @param order Order of the rational approximation
#' @param nu Smoothness parameter
#' @param kappa Range parameter
#' @param sigma Standard deviation
#' @param type_rational Method used to compute the coefficients of the rational approximation.
#' @param type_interp Interpolation method for the rational coefficients. 
#' @param equally_spaced Are the locations equally spaced?
#' @param ordering Return the matrices ordered by field or by location? 
#' @return A list containing the precision matrix `Q` of the process and its derivatives if they exist, and
#' a matrix `A` that extracts the elements corresponding to the process. 
#' @export 
#'
#' @examples
#' h <- seq(from = 0, to = 1, length.out = 100)
#' cov.true <- matern.covariance(h, kappa = 10, sigma = 1, nu = 0.8)
#' Q.approx <- matern.rational.precision(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)
#' cov.approx <- Q.approx$A%*%solve(Q.approx$Q,Q.approx$A[1,])
#' 
#' plot(h, cov.true)
#' lines(h, cov.approx, col = 2)
matern.rational.precision <- function(loc,
                                      order,
                                      nu,
                                      kappa,
                                      sigma,
                                      type_rational = "brasil",
                                      type_interp = "spline",
                                      equally_spaced = FALSE,
                                      ordering = c("field", "location")) {
    ordering <- ordering[[1]]
    if(!(ordering %in% c("field", "location"))) {
        stop("Ordering must be 'field' or 'location'.")
    }
    if(is.matrix(loc) && min(dim(h)) > 1) {
        stop("Only one dimensional locations supported.")
    }
    alpha=nu+1/2
    n <- length(loc)
    if (nu >= 3/2) {
        stop("Only implemented for 0 < nu < 3/2")
    } else if(nu==1/2){
        scaling <- gamma(nu)/(sigma^2*gamma(alpha)*sqrt(4*pi)*kappa^(2*nu))
        Q <- exp_precision(loc,kappa)
        A <- Diagonal(n)
    } else {
        coeff <- interp_rational_coefficients(order = order, 
                                              type_rational_approx = type_rational, 
                                              type_interp = type_interp, 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
        
        A <- I <- Diagonal(n)
        scaling <- kappa^(2*alpha)*gamma(nu)/(sigma^2*gamma(alpha)*sqrt(4*pi)*kappa^(2*nu))
        
        if(nu < 1/2) {
            Q <- I/k
            for(i in 1:length(p)){
                Q <- bdiag(Q, exp_precision(loc,sqrt(kappa^2*(1-p[i])))/(r[i]*kappa^2))
                A <- cbind(A, I)
            }
        } else {
            Q <- exp_precision(loc,kappa)/(k*kappa^2)
            
            for(i in 1:length(p)){
                Q_tmp <- matern.p.precision(loc, kappa, kappa^2*p[i],
                                            equally_spaced = equally_spaced, alpha)/(r[i]*kappa^4)
                Q <- bdiag(Q, Q_tmp)
                A = cbind(A,kronecker(Diagonal(n),Matrix(c(1,0),nrow=1,ncol=2)))
            }
        }
    }
    if(ordering == "location") {
        reo <- compute.reordering(n,order,alpha)
        Q <- Q[reo,reo]
        A <- A[,reo]
    } 
    return(list(Q = scaling*Q, A = A))    
}

### Utility functions below ####

# Derivatives of the matern covariance
matern.derivative = function(h, kappa, nu, sigma, deriv = 1) 
{
    if (deriv == 1) {
        C = h * matern.covariance(h, kappa = kappa, nu = nu - 
                                      1, sigma = sigma)
        C[h == 0] = 0
    }
    else if (deriv == 2) {
        C = matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma) + 
            h * matern_derivative(h, kappa = kappa, nu = nu - 1, sigma = sigma, deriv = 1)
    }
    else {
        C = (deriv - 1) * matern_derivative(h, kappa = kappa, 
                                            nu = nu - 1, sigma = sigma, deriv = deriv - 2) + 
            h * matern_derivative(h, kappa = kappa, nu = nu - 
                                      1, sigma = sigma, deriv = deriv - 1)
    }
    return(-(kappa^2/(2 * (nu - 1))) * as.matrix(C))
}



# Reorder matern
compute.reordering <- function(n,m,alpha) {
    if(alpha <0 || alpha > 2) {
        stop("only 0<alpha<2 supported")
    }
    if(alpha < 1) {
        return(as.vector(rep(seq(from=1,to=(m+1)*n,by=n),n) + kronecker(0:(n-1),rep(1,m+1))))
    } else {
        tmp <- rep(c(1,c(rbind(1+n*seq(from=1,to=2*m,by=2), 2+n*seq(from=1,to=2*m,by=2)))),n) 
        return(tmp + as.vector(rbind(0:(n-1), matrix(rep(1,2*m*n),2*m,n)%*%Diagonal(n,2*c(0:(n-1))))))
    }
}

# Precision for exponential covariance
exp_precision <- function(loc, kappa, boundary = "free") {
    n <- length(loc)
    l <- diff(loc)
    
    # Precompute reusable terms
    a <- exp(-2 * kappa * l)
    b <- 1 - a
    
    # Initialize matrix Q efficiently
    Q <- Matrix(0, nrow = n, ncol = n)
    
    # Set diagonal elements
    diag(Q)[1:(n-1)] <- 1/2 + a/b
    diag(Q)[2:n] <- diag(Q)[2:n] + 1/2 + a/b
    
    # Set off-diagonal elements
    #  indices <- c(which(row(Q) == col(Q) + 1), which(row(Q) == col(Q) - 1))
    # off_diag_values <- -exp(-kappa * l) /b
    # Q[indices] <- off_diag_values
    
    
    row_indices <- seq_len(n - 1)
    col_indices <- row_indices + 1  # Offsets for off-diagonal elements
    
    # Compute the off-diagonal values
    off_diag_values <- -exp(-kappa * l)/b
    
    # Set the off-diagonal elements in matrix Q
    Q[cbind(row_indices, col_indices)] <- off_diag_values
    Q[cbind(col_indices, row_indices)] <- off_diag_values
    
    # Adjust boundary conditions if required
    if (boundary == "free") {
        Q[1, 1] <- Q[1, 1] + 0.5
        Q[n, n] <- Q[n, n] + 0.5
    }
    
    return(2 * kappa * Q[])
}


matern.p.precision <- function(loc,kappa,p,equally_spaced = FALSE, alpha = 1){
    fa <- floor(alpha)
    n <- length(loc)
    #Q <- Matrix(0, nrow = (fa + 1) * n, ncol = (fa + 1) * n, sparse = TRUE)
    # number of non-zero entries: 12(n-2) + 16 for n >= 2.
    # N <- 12*(n-2) + 16
    #n*4 + (n-1)*4 = 8n - 4 - n    
    #n*((fa+1)^2) + (n-1)*(fa+1)^2 - 1 6 (n-1)+(n-2)+...+1 = (n-1)*n/2, (2fa -1)*fa
    #N <- (4*fa^2)*n + (4*fa^2)*(n-1)- n*(2fa -1)*fa
    N <- 6*n*fa^2 + n*fa - 4*fa^2 
    #N <- 7*n-4
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    
    if(equally_spaced){
        Q.1 <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p), matern.p.joint(loc[1],loc[1+1],kappa,p,fa)),
                           cbind(matern.p.joint(loc[1+1],loc[1],kappa,p), matern.p.joint(loc[1+1],loc[1+1],kappa,p,fa))))
        
        Q.i <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,fa), 
                                 matern.p.joint(loc[1],loc[2],kappa,p,fa),
                                 matern.p.joint(loc[1],loc[3],kappa,p,fa)),
                           cbind(matern.p.joint(loc[2],loc[1],kappa,p,fa),
                                 matern.p.joint(loc[2],loc[2],kappa,p,fa),
                                 matern.p.joint(loc[2],loc[3],kappa,p,fa)),
                           cbind(matern.p.joint(loc[3],loc[1],kappa,p,fa),
                                 matern.p.joint(loc[3],loc[2],kappa,p,fa),
                                 matern.p.joint(loc[3],loc[3],kappa,p,fa))))[-c(1:(2*fa)),-c(1:(2*fa))]
    }
    
    
    for(i in 1:(n-1)){
        if(i==1){
            if(!equally_spaced){
                Q.i <- solve(rbind(cbind(matern.p.joint(loc[i],loc[i],kappa,p,fa),
                                         matern.p.joint(loc[i],loc[i+1],kappa,p,fa)),
                                   cbind(matern.p.joint(loc[i+1],loc[i],kappa,p,fa),
                                         matern.p.joint(loc[i+1],loc[i+1],kappa,p,fa))))
                counter <- 1
                for(ki in 1:(2*fa)){
                    for(kj in ki:(4*fa)){
                        ii[counter] <- ki
                        jj[counter] <- kj
                        val[counter] <- Q.i[ki,kj]
                        counter <- counter + 1
                    }
                }
            } else{
                counter <- 1
                for(ki in 1:(2*fa)) {
                    for(kj in ki:(4*fa)) {
                        ii[counter] <- ki
                        jj[counter] <- kj
                        val[counter] <- Q.1[ki,kj]
                        counter <- counter + 1
                    }
                }
            }
            
        } else {
            if(!equally_spaced){
                Q.i <- solve(rbind(cbind(matern.p.joint(loc[i-1],loc[i-1],kappa,p,fa),
                                         matern.p.joint(loc[i-1],loc[i],kappa,p,fa),
                                         matern.p.joint(loc[i-1],loc[i+1],kappa,p,fa)),
                                   cbind(matern.p.joint(loc[i],loc[i-1],kappa,p,fa),
                                         matern.p.joint(loc[i],loc[i],kappa,p,fa),
                                         matern.p.joint(loc[i],loc[i+1],kappa,p,fa)),
                                   cbind(matern.p.joint(loc[i+1],loc[i-1],kappa,p,fa),
                                         matern.p.joint(loc[i+1],loc[i],kappa,p,fa),
                                         matern.p.joint(loc[i+1],loc[i+1],kappa,p,fa))))[-c(1:(2*fa)),-c(1:(2*fa))]
            } 
            
            for(ki in 1:(2*fa)){
                for(kj in ki:(4*fa)){
                    ii[counter] <- 2*fa*i-1 + ki-1
                    jj[counter] <- 2*fa*i-1 + kj-1
                    val[counter] <-Q.i[ki,kj]
                    counter <- counter + 1
                }
            }         
        }
    }
    for(ki in 1:(2*fa)){
        for(kj in (2+ki):(4*fa)){
            ii[counter] <- 2*fa*n-1 + (ki-1)
            jj[counter] <- 2*fa*n-3 + kj-1
            val[counter] <-Q.i[ki+2,kj]
            counter <- counter + 1
        }
    }    

    Q <- Matrix::sparseMatrix(i   = ii,
                              j    = jj,
                              x    = val,
                              dims = c(2*fa*n, 2*fa*n), symmetric=TRUE)
    return(Q)
}


#Joint covariance of process and derivative for shifted Matern
matern.p.joint <- function(s,t,kappa,p, fa = 1){
    mat <- matrix(0, nrow = 2*fa, ncol = 2*fa)
    for(i in 1:(2*fa)) {
        for(j in i:(2*fa)) {
            if(i==j) {
                mat[i,j] <- ((-1)^(i-1))*matern.p(s,t,kappa,p,fa, deriv = 2*(i-1))
            } else {
                mat[i,j] <- matern.p(s,t,kappa,p,fa, deriv = i-1 + j - 1)
                mat[j,i] <- -mat[i,j]    
            }
            
        }
    }
    return(mat)
}


#covariance for shifted Matern for 1<alpha<2
matern.p <- function(s,t,kappa,p,fa, deriv=0){
    h <- s-t
    # print(p)
    if(p==0){
        if(deriv==0){
            return(matern.covariance(h, kappa = kappa, nu = fa + 1/2, sigma = 1))
        } else {
            return(matern.derivative(h, kappa = kappa, nu = fa + 1/2, sigma = 1, deriv = deriv))
        }
    } else {
        if(fa == 1) {
            a <- -1/(2*kappa*p)
            b <- sqrt(kappa^2-p)
            if(deriv==0){
                return(a*(exp(-kappa*abs(h))-(kappa/b)*exp(-b*abs(h))))
            } else if(deriv == 1) {
                return(-a*kappa*sign(h)*(-exp(-kappa*abs(h)) + exp(-b*abs(h))))
            } else if(deriv == 2) {
                return(a*kappa*(kappa*exp(-kappa*abs(h))-b*exp(-b*abs(h))))
            } else{
                stop("only deriv=0,1,2 allowed")
            }    
        } else {
          stop("not implemented yet")   
        }
    }
}


