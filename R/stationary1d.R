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
    
    if(is.matrix(h) && min(dim(h)) > 1) {
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
    if (alpha%%1 == 0) {
        tmp <- matern.p.precision(loc = loc,kappa = kappa,p = 0,
                                  equally_spaced = equally_spaced, alpha = alpha) 
        Q <- tmp$Q
        A <- tmp$A
    } else {
        coeff <- interp_rational_coefficients(order = order, 
                                              type_rational_approx = type_rational, 
                                              type_interp = type_interp, 
                                              alpha = alpha)
        r <- coeff$r
        p <- coeff$p
        k <- coeff$k
        
        ## k part
        tmp <- matern.k.precision(loc = loc,kappa,equally_spaced = equally_spaced, 
                                  alpha = alpha)
        Q <- tmp$Q/(k*sigma^2)
        A <- tmp$A
        
        ## p part
        for(i in 1:length(p)){
            tmp <- matern.p.precision(loc = loc, kappa = kappa, p =p[i],
                                      equally_spaced = equally_spaced, 
                                      alpha = alpha)
            Q <- bdiag(Q, tmp$Q/(r[i]*sigma^2))
            A = cbind(A,tmp$A)
        }
    }
    if(ordering == "location") {
        reo <- compute.reordering(n,order,alpha)
        Q <- Q[reo,reo]
        A <- A[,reo]
    } 
    return(list(Q = Q, A = A))    
}

### Utility below

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

# Derivatives of the matern covariance
matern.derivative = function(h, kappa, nu, sigma, deriv = 1) 
{
    if(deriv == 0) {
        C = matern.covariance(h, kappa = kappa, nu = nu, sigma = sigma)
    } else if (deriv == 1) {
        C = h * matern.covariance(h, kappa = kappa, nu = nu - 
                                      1, sigma = sigma)
        C[h == 0] = 0
    }
    else if (deriv == 2) {
        C = matern.covariance(h, kappa = kappa, nu = nu - 1, sigma = sigma) + 
            h * matern.derivative(h, kappa = kappa, nu = nu - 1, sigma = sigma, deriv = 1)
    }
    else {
        C = (deriv - 1) * matern.derivative(h, kappa = kappa, 
                                            nu = nu - 1, sigma = sigma, deriv = deriv - 2) + 
            h * matern.derivative(h, kappa = kappa, nu = nu - 
                                      1, sigma = sigma, deriv = deriv - 1)
    }
    return(-(kappa^2/(2 * (nu - 1))) * as.matrix(C))
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







#Joint covariance of process and derivative for shifted Matern
matern.p.joint <- function(s,t,kappa,p, alpha = 1){
    
    if(alpha%%1 == 0) {
        fa <- alpha
    } else {
        fa <- floor(alpha) + 1    
    }
    mat <- matrix(0, nrow = fa, ncol = fa)
    for(i in 1:fa) {
        for(j in i:fa) {
            if(i==j) {
                mat[i,i] <- ((-1)^(i-1))*matern.p.deriv(s,t,kappa,p,alpha, deriv = 2*(i-1))
            } else {
                mat[i,j] <- (-1)^(j-1)*matern.p.deriv(s,t,kappa,p,alpha, deriv = i-1 + j - 1)
                mat[j,i] <- (-1)^(i-1)*matern.p.deriv(s,t,kappa,p,alpha, deriv = i-1 + j - 1)
            }
        }
    }
    
    return(mat)
}


matern.p <- function(s,t,kappa,p,alpha){
    h <- s-t
    if(p==0){
        return(matern.covariance(h, kappa = kappa, nu = alpha - 1/2, sigma = 1))
    } else {
        ca <- gamma(alpha)/gamma(alpha-0.5)
        fa <- floor(alpha)
        kp <- kappa*sqrt(1-p)
        out <- matern.covariance(h, kappa = kp, nu = 1/2, 
                                 sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)))
        if(alpha < 1) {
            return(out)
        } else {
            out <- out/p^fa
            for(j in 1:fa) {
                out <- out - matern.covariance(h, kappa = kappa, nu = j-1/2, 
                                               sigma = sqrt(ca*gamma(j-0.5)/gamma(j)))/p^(fa + 1 - j)
            }
            return(out)    
        }
    }
}

matern.p.deriv <- function(s,t,kappa,p,alpha,deriv = 0){
    h <- s-t
    if(deriv ==0){
        return(matern.p(s,t,kappa,p,alpha))
    } else {
        if(p==0){
            return(matern.derivative(h, kappa = kappa, nu = alpha - 1/2, 
                                     sigma = 1, deriv = deriv))
        } else {
            ca <- gamma(alpha)/gamma(alpha-0.5)
            fa <- floor(alpha)
            kp <- kappa*sqrt(1-p)
            out <- matern.derivative(h, kappa = kp, nu =  1/2, 
                                     sigma = sqrt(ca*sqrt(pi)/sqrt(1-p)), deriv = deriv)
            if(alpha < 1) {
                return(out)
            } else {
                out <- out/p^fa
                for(j in 1:fa) {
                    out <- out - matern.derivative(h, kappa = kappa, nu = j-1/2, 
                                                   sigma = sqrt(ca*gamma(j-0.5)/gamma(j)), 
                                                   deriv = deriv)/p^(fa + 1 - j)
                }
            }
            return(out)
        }    
    }
    
}

matern.p.precision <- function(loc,kappa,p,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    if(alpha==1) {
        Q <- exp_precision(loc,kappa)/(2*kappa)
        A <- Diagonal(n)
        
        return(list(Q=Q,A=A))
    } else {
        if(alpha%%1 == 0) {
            fa <- alpha
        } else {
            fa <- floor(alpha) + 1    
        }
        
        if(fa == 1) {
            N <- n  + n - 1 
        } else {
            N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
        }
        
        ii <- numeric(N)
        jj <- numeric(N)
        val <- numeric(N)
        
        if(equally_spaced){
            Q.1 <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1],loc[1+1],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[1+1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1+1],loc[1+1],kappa,p,alpha))))
            
            Q.i <- solve(rbind(cbind(matern.p.joint(loc[1],loc[1],kappa,p,alpha), 
                                     matern.p.joint(loc[1],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[1],loc[3],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[2],loc[1],kappa,p,alpha),
                                     matern.p.joint(loc[2],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[2],loc[3],kappa,p,alpha)),
                               cbind(matern.p.joint(loc[3],loc[1],kappa,p,alpha),
                                     matern.p.joint(loc[3],loc[2],kappa,p,alpha),
                                     matern.p.joint(loc[3],loc[3],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
        }
        
        
        for(i in 1:(n-1)){
            if(i==1){
                if(!equally_spaced){
                    Q.1 <- solve(rbind(cbind(matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))
                    
                } 
                counter <- 1
                for(ki in 1:fa) {
                    for(kj in ki:(2*fa)) {
                        ii[counter] <- ki
                        jj[counter] <- kj
                        val[counter] <- Q.1[ki,kj]
                        counter <- counter + 1
                    }
                }
            } else {
                if(!equally_spaced){
                    Q.i <- solve(rbind(cbind(matern.p.joint(loc[i-1],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i-1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i-1],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i],loc[i+1],kappa,p,alpha)),
                                       cbind(matern.p.joint(loc[i+1],loc[i-1],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i],kappa,p,alpha),
                                             matern.p.joint(loc[i+1],loc[i+1],kappa,p,alpha))))[-c(1:fa),-c(1:fa)]
                } 
                # Q[(2*n-1):(2*n), (2*n-3):(2*n)] = Q.i[3:4,]
                for(ki in 1:fa){
                    for(kj in ki:(2*fa)){
                        ii[counter] <- fa*(i-1) + ki
                        jj[counter] <- fa*(i-1) + kj
                        val[counter] <-Q.i[ki,kj]
                        counter <- counter + 1
                    }
                }     
            }
        }
        for(ki in 1:fa){
            for(kj in ki:fa){
                ii[counter] <- fa*(n-1) + ki
                jj[counter] <- fa*(n-1) + kj
                val[counter] <-Q.i[ki+fa,kj+fa]
                counter <- counter + 1
            }
        }    
        Q <- Matrix::sparseMatrix(i   = ii,
                                  j    = jj,
                                  x    = val,
                                  dims = c(fa*n, fa*n), symmetric=TRUE)
        
        A <-  Matrix::sparseMatrix(i   = 1:n,
                                   j    = seq(from=1,to=n*fa,by=fa),
                                   x    = rep(1,n),
                                   dims = c(n, fa*n))
        return(list(Q=Q,A=A))    
    }
}

matern.k.precision <- function(loc,kappa,equally_spaced = FALSE, alpha = 1) {
    
    n <- length(loc)
    if(alpha==1) {
        Q <- exp_precision(loc,kappa)/(2*kappa)
        A <- Diagonal(n)
        return(list(Q=Q,A=A))
    } else {
        fa <- floor(alpha)
        if(fa == 0) {
            N <- n 
        } else {
            N <- n*fa^2 + (n-1)*fa^2 - n*fa*(fa -1)/2    
        }
        da <- max(fa,1)
        ii <- numeric(N)
        jj <- numeric(N)
        val <- numeric(N)
        
        if(equally_spaced){
            Q.1 <- solve(rbind(cbind(matern.k.joint(loc[1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1],loc[1+1],kappa,alpha)),
                               cbind(matern.k.joint(loc[1+1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1+1],loc[1+1],kappa,alpha))))
            
            Q.i <- solve(rbind(cbind(matern.k.joint(loc[1],loc[1],kappa,alpha), 
                                     matern.k.joint(loc[1],loc[2],kappa,alpha),
                                     matern.k.joint(loc[1],loc[3],kappa,alpha)),
                               cbind(matern.k.joint(loc[2],loc[1],kappa,alpha),
                                     matern.k.joint(loc[2],loc[2],kappa,alpha),
                                     matern.k.joint(loc[2],loc[3],kappa,alpha)),
                               cbind(matern.k.joint(loc[3],loc[1],kappa,alpha),
                                     matern.k.joint(loc[3],loc[2],kappa,alpha),
                                     matern.k.joint(loc[3],loc[3],kappa,alpha))))[-c(1:da),-c(1:da)]
            
        }
    }
    
    for(i in 1:(n-1)){
        if(i==1){
            if(!equally_spaced){
                Q.1 <- solve(rbind(cbind(matern.k.joint(loc[i],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i+1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i+1],kappa,alpha))))
                
            } 
            counter <- 1
            for(ki in 1:da) {
                for(kj in ki:(2*da)) {
                    ii[counter] <- ki
                    jj[counter] <- kj
                    val[counter] <- Q.1[ki,kj]
                    counter <- counter + 1
                }
            }
        } else {
            if(!equally_spaced){
                Q.i <- solve(rbind(cbind(matern.k.joint(loc[i-1],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i-1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i-1],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i],loc[i+1],kappa,alpha)),
                                   cbind(matern.k.joint(loc[i+1],loc[i-1],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i],kappa,alpha),
                                         matern.k.joint(loc[i+1],loc[i+1],kappa,alpha))))[-c(1:da),-c(1:da)]
                
            } 
            # Q[(2*n-1):(2*n), (2*n-3):(2*n)] = Q.i[3:4,]
            for(ki in 1:da){
                for(kj in ki:(2*da)){
                    ii[counter] <- da*(i-1) + ki
                    jj[counter] <- da*(i-1) + kj
                    val[counter] <-Q.i[ki,kj]
                    counter <- counter + 1
                }
            }     
        }
    }
    for(ki in 1:da){
        for(kj in ki:da){
            ii[counter] <- da*(n-1) + ki
            jj[counter] <- da*(n-1) + kj
            val[counter] <-Q.i[ki+da,kj+da]
            counter <- counter + 1
        }
    }    
    Q <- Matrix::sparseMatrix(i   = ii,
                              j    = jj,
                              x    = val,
                              dims = c(da*n, da*n), symmetric=TRUE)
    
    A <-  Matrix::sparseMatrix(i   = 1:n,
                               j    = seq(from=1,to=n*da,by=da),
                               x    = rep(1,n),
                               dims = c(n, da*n))
    return(list(Q=Q,A=A))    
}

#Joint covariance of process and derivative for shifted Matern
matern.k.joint <- function(s,t,kappa,alpha = 1){
    fa <- floor(alpha)
    if(fa == 0) {
        mat <- matrix(0, nrow = 1, ncol = 1)
        mat[1,1] <- matern.k(s,t,kappa,alpha)
    } else {
        mat <- matrix(0, nrow = fa, ncol = fa)    
        for(i in 1:fa) {
            for(j in i:fa) {
                if(i==j) {
                    mat[i,i] <- ((-1)^(i-1))*matern.k.deriv(s,t,kappa,alpha, deriv = 2*(i-1))
                } else {
                    mat[i,j] <- (-1)^(j-1)*matern.k.deriv(s,t,kappa,alpha, deriv = i-1 + j - 1)
                    mat[j,i] <- (-1)^(i-1)*matern.k.deriv(s,t,kappa,alpha, deriv = i-1 + j - 1)
                }
            }
        }
    }
    return(mat)
}

matern.k.deriv <- function(s,t,kappa,alpha, deriv = 0){
    h <- s-t
    if(alpha<1){
        if(deriv == 0) {
            return((h==0)*sqrt(4*pi)*ca/kappa)    
        } else {
            stop("Not differentiable")
        }
    } else {
        fa <- floor(alpha)
        ca <- gamma(alpha)/gamma(alpha-0.5)
        cfa <- gamma(fa)/gamma(fa-0.5)
        return(matern.derivative(h, kappa = kappa, nu = fa - 1/2, 
                                 sigma = sqrt(ca/cfa), deriv = deriv))
    }
}


matern.k <- function(s,t,kappa,alpha){
    h <- s-t
    
    if(alpha<1){
        return((h==0)*sqrt(4*pi)*ca/kappa)
    } else {
        fa <- floor(alpha)
        ca <- gamma(alpha)/gamma(alpha-0.5)
        cfa <- gamma(fa)/gamma(fa-0.5)
        return(matern.covariance(h, kappa = kappa, nu = fa - 1/2, sigma = sqrt(ca/cfa)))
    }
}
