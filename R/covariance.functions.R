#' The Matern covariance function
#'
#' `matern.covariance` evaluates the Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h).}
#'
#' @param h Distances to evaluate the covariance function at.
#' @param kappa Range parameter.
#' @param nu Shape parameter.
#' @param sigma Standard deviation.
#'
#' @return A vector with the values C(h).
#' @export
#'
#' @examples
#' x <- seq(from = 0, to = 1, length.out = 101)
#' plot(x, matern.covariance(abs(x - 0.5), kappa = 10, nu = 1 / 5, sigma = 1),
#'   type = "l", ylab = "C(h)", xlab = "h"
#' )
#'
matern.covariance <- function(h,
                              kappa,
                              nu,
                              sigma) {
    if (nu == 1 / 2) {
        C <- sigma^2 * exp(-kappa * abs(h))
    } else {
        C <- (sigma^2 / (2^(nu - 1) * gamma(nu))) *
            ((kappa * abs(h))^nu) * besselK(kappa * abs(h), nu)
        C[h == 0] <- sigma^2
    }
    
    #return(as.matrix(C))
    return(C)
}

#' The 1d folded Matern covariance function
#'
#' @description
#' `folded.matern.covariance.1d` evaluates the 1d
#' folded Matern covariance function over an interval \eqn{[0,L]}.
#'
#' @details
#' `folded.matern.covariance.1d` evaluates the 1d folded Matern
#' covariance function over an interval \eqn{[0,L]} under different
#' boundary conditions. For periodic boundary conditions
#' \deqn{C_{\mathcal{P}}(h,m) = \sum_{k=-\infty}^{\infty} (C(h-m+2kL),}
#' for Neumann boundary conditions
#' \deqn{C_{\mathcal{N}}(h,m) = \sum_{k=-\infty}^{\infty}
#' (C(h-m+2kL)+C(h+m+2kL)),}
#' and for Dirichlet boundary conditions:
#' \deqn{C_{\mathcal{D}}(h,m) = \sum_{k=-\infty}^{\infty}
#' (C(h-m+2kL)-C(h+m+2kL)),}
#' where \eqn{C(\cdot)} is the Matern covariance function:
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu
#' K_\nu(\kappa h).}
#'
#' We consider the truncation:
#' \deqn{C_{{\mathcal{P}},N}(h,m) = \sum_{k=-N}^{N} C(h-m+2kL),
#' C_{\mathcal{N},N}(h,m) = \sum_{k=-\infty}^{\infty}
#' (C(h-m+2kL)+C(h+m+2kL)),}
#' and
#' \deqn{C_{\mathcal{D},N}(h,m) = \sum_{k=-N}^{N}
#' (C(h-m+2kL)-C(h+m+2kL)).}
#'
#' @param h,m Vectors of arguments of the covariance function.
#' @param kappa Range parameter.
#' @param nu Shape parameter.
#' @param sigma Standard deviation.
#' @param L The upper bound of the interval \eqn{[0,L]}. By default, `L=1`.
#' @param N The truncation parameter.
#' @param boundary The boundary condition. The possible conditions
#' are `"neumann"` (default), `"dirichlet"` or
#' `"periodic"`.
#'
#' @return A matrix with the corresponding covariance values.
#' @export
#'
#' @examples
#' x <- seq(from = 0, to = 1, length.out = 101)
#' plot(x, folded.matern.covariance.1d(rep(0.5, length(x)), x,
#'   kappa = 10, nu = 1 / 5, sigma = 1
#' ),
#' type = "l", ylab = "C(h)", xlab = "h"
#' )
#'
folded.matern.covariance.1d <- function(h, m, kappa, nu, sigma,
                                        L = 1, N = 10,
                                        boundary = c(
                                            "neumann",
                                            "dirichlet", "periodic"
                                        )) {
    boundary <- tolower(boundary[1])
    if (!(boundary %in% c("neumann", "dirichlet", "periodic"))) {
        stop("The possible boundary conditions are 'neumann',
    'dirichlet' or 'periodic'!")
    }
    if (length(h) != length(m)) {
        stop("h and m should have the same length!")
    }
    
    s1 <- sapply(-N:N, function(j) {
        h - m + 2 * j * L
    })
    s2 <- sapply(-N:N, function(j) {
        h + m + 2 * j * L
    })
    if (boundary == "neumann") {
        C <- rowSums(matern.covariance(
            h = s1, kappa = kappa,
            nu = nu, sigma = sigma
        ) +
            matern.covariance(
                h = s2, kappa = kappa,
                nu = nu, sigma = sigma
            ))
    } else if (boundary == "dirichlet") {
        C <- rowSums(matern.covariance(
            h = s1, kappa = kappa,
            nu = nu, sigma = sigma
        ) -
            matern.covariance(
                h = s2, kappa = kappa,
                nu = nu, sigma = sigma
            ))
    } else {
        C <- rowSums(matern.covariance(
            h = s1,
            kappa = kappa, nu = nu, sigma = sigma
        ))
    }
    
    if (length(h) == 1) {
        return(sum(C))
    }
    return(as.matrix(C))
}

#' The 2d folded Matern covariance function
#'
#' @description
#' `folded.matern.covariance.2d` evaluates the 2d
#' folded Matern covariance function over an interval
#' \eqn{[0,L]\times [0,L]}.
#'
#' @details
#' `folded.matern.covariance.2d` evaluates the 1d folded
#' Matern covariance function over an interval
#' \eqn{[0,L]\times [0,L]} under different boundary conditions.
#' For periodic boundary conditions
#' \deqn{C_{\mathcal{P}}((h_1,h_2),(m_1,m_2)) =
#' \sum_{k_2=-\infty}^\infty \sum_{k_1=-\infty}^{\infty}
#' (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|),}
#' for Neumann boundary conditions
#' \deqn{C_{\mathcal{N}}((h_1,h_2),(m_1,m_2)) =
#' \sum_{k_2=-\infty}^\infty \sum_{k_1=-\infty}^{\infty}
#' (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|)+C(\|(h_1-m_1+2k_1L,
#' h_2+m_2+2k_2L)\|)+C(\|(h_1+m_1+2k_1L,h_2-m_2+2k_2L)\|)+
#' C(\|(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\|)),}
#' and for Dirichlet boundary conditions:
#' \deqn{C_{\mathcal{D}}((h_1,h_2),(m_1,m_2)) = \sum_{k_2=-\infty}^\infty
#' \sum_{k_1=-\infty}^{\infty} (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|)-
#' C(\|(h_1-m_1+2k_1L,h_2+m_2+2k_2L)\|)-C(\|(h_1+m_1+2k_1L,
#' h_2-m_2+2k_2L)\|)+C(\|(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\|)),}
#' where \eqn{C(\cdot)} is the Matern covariance function:
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}
#' (\kappa h)^\nu K_\nu(\kappa h).}
#'
#' We consider the truncation for \eqn{k_1,k_2} from \eqn{-N} to \eqn{N}.
#'
#' @param h,m Vectors with two coordinates.
#' @param kappa Range parameter.
#' @param nu Shape parameter.
#' @param sigma Standard deviation.
#' @param L The upper bound of the square \eqn{[0,L]\times [0,L]}.
#' By default, `L=1`.
#' @param N The truncation parameter.
#' @param boundary The boundary condition. The possible conditions
#' are `"neumann"` (default), `"dirichlet"`,
#' `"periodic"` or `"R2"`.
#'
#' @return The correspoding covariance.
#' @export
#'
#' @examples
#' h <- c(0.5, 0.5)
#' m <- c(0.5, 0.5)
#' folded.matern.covariance.2d(h, m, kappa = 10, nu = 1 / 5, sigma = 1)
#'
folded.matern.covariance.2d <- function(h, m, kappa, nu, sigma,
                                        L = 1, N = 10,
                                        boundary = c(
                                            "neumann",
                                            "dirichlet",
                                            "periodic",
                                            "R2"
                                        )) {
    boundary <- tolower(boundary[1])
    if (!(boundary %in% c(
        "neumann", "dirichlet",
        "periodic", "r2"
    ))) {
        stop("The possible boundary conditions are
    'neumann', 'dirichlet', 'periodic' or 'R2'!")
    }
    
    if (is.vector(h)) {
        if (!is.vector(m)) {
            stop("since 'h' is a vector, 'm' should be a vector!")
        }
        
        if ((length(h) != 2) || (length(m) != 2)) {
            stop("The vectors h and m should have length 2!")
        }
    } else if (is.matrix(h)) {
        if (!is.matrix(m)) {
            stop("since 'h' is a matrix, 'm' should be a matrix!")
        }
        if (ncol(h) != 2) {
            stop("h must have two columns!")
        }
        if (!all(dim(h) == dim(m))) {
            stop("h and m must have the same dimensions!")
        }
    }
    
    if (!is.vector(h) && !is.matrix(h)) {
        stop("h should be either a vector or a matrix!")
    }
    
    list.comb <- expand.grid(-N:N, -N:N)
    
    if (is.matrix(h)) {
        h_matrix_1 <- matrix(rep(h[, 1], length(list.comb[, 1])), nrow = nrow(h))
        h_matrix_2 <- matrix(rep(h[, 2], length(list.comb[, 1])), nrow = nrow(h))
        m_matrix_1 <- matrix(rep(m[, 1], length(list.comb[, 1])), nrow = nrow(m))
        m_matrix_2 <- matrix(rep(m[, 2], length(list.comb[, 1])), nrow = nrow(m))
        list_comb_matrix_1 <- t(matrix(rep(
            list.comb[, 1],
            nrow(h)
        ), ncol = nrow(h)))
        list_comb_matrix_2 <- t(matrix(rep(
            list.comb[, 2],
            nrow(h)
        ), ncol = nrow(h)))
        
        s11 <- sqrt((h_matrix_1 - m_matrix_1 + 2 *
                         list_comb_matrix_1 * L)^2 + (h_matrix_2 -
                                                          m_matrix_2 + 2 * list_comb_matrix_2 * L)^2)
        s12 <- sqrt((h_matrix_1 - m_matrix_1 + 2 *
                         list_comb_matrix_1 * L)^2 + (h_matrix_2 +
                                                          m_matrix_2 + 2 * list_comb_matrix_2 * L)^2)
        s21 <- sqrt((h_matrix_1 + m_matrix_1 + 2 *
                         list_comb_matrix_1 * L)^2 + (h_matrix_2 -
                                                          m_matrix_2 + 2 * list_comb_matrix_2 * L)^2)
        s22 <- sqrt((h_matrix_1 + m_matrix_1 + 2 *
                         list_comb_matrix_1 * L)^2 + (h_matrix_2 +
                                                          m_matrix_2 + 2 * list_comb_matrix_2 * L)^2)
        
        if (boundary == "neumann") {
            C <- rowSums(matern.covariance(
                h = s11,
                kappa = kappa, nu = nu, sigma = sigma
            ) +
                matern.covariance(
                    h = s12, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s21, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s22, kappa = kappa,
                    nu = nu, sigma = sigma
                ))
        } else if (boundary == "dirichlet") {
            C <- rowSums(matern.covariance(
                h = s11,
                kappa = kappa, nu = nu, sigma = sigma
            ) -
                matern.covariance(
                    h = s12, kappa = kappa,
                    nu = nu, sigma = sigma
                ) -
                matern.covariance(
                    h = s21, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s22, kappa = kappa,
                    nu = nu, sigma = sigma
                ))
        } else if (boundary == "r2") {
            C <- matern.covariance(
                h = sqrt((h[1] - m[1])^2 +
                             (h[2] - m[2])^2), kappa = kappa, sigma = sigma,
                nu = nu
            )
        } else {
            C <- rowSums(matern.covariance(
                h = s11,
                kappa = kappa, nu = nu, sigma = sigma
            ))
        }
        return(as.double(C))
    } else {
        s11 <- sqrt((h[1] - m[1] + 2 * list.comb[, 1] * L)^2 +
                        (h[2] - m[2] + 2 * list.comb[, 2] * L)^2)
        s12 <- sqrt((h[1] - m[1] + 2 * list.comb[, 1] * L)^2 +
                        (h[2] + m[2] + 2 * list.comb[, 2] * L)^2)
        s21 <- sqrt((h[1] + m[1] + 2 * list.comb[, 1] * L)^2 +
                        (h[2] - m[2] + 2 * list.comb[, 2] * L)^2)
        s22 <- sqrt((h[1] + m[1] + 2 * list.comb[, 1] * L)^2 +
                        (h[2] + m[2] + 2 * list.comb[, 2] * L)^2)
        
        if (boundary == "neumann") {
            C <- sum(matern.covariance(
                h = s11, kappa = kappa,
                nu = nu, sigma = sigma
            ) +
                matern.covariance(
                    h = s12, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s21, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s22, kappa = kappa,
                    nu = nu, sigma = sigma
                ))
        } else if (boundary == "dirichlet") {
            C <- sum(matern.covariance(
                h = s11,
                kappa = kappa, nu = nu, sigma = sigma
            ) -
                matern.covariance(
                    h = s12, kappa = kappa,
                    nu = nu, sigma = sigma
                ) -
                matern.covariance(
                    h = s21, kappa = kappa,
                    nu = nu, sigma = sigma
                ) +
                matern.covariance(
                    h = s22, kappa = kappa,
                    nu = nu, sigma = sigma
                ))
        } else if (boundary == "r2") {
            C <- matern.covariance(
                h = sqrt((h[1] - m[1])^2 +
                             (h[2] - m[2])^2), kappa = kappa, sigma = sigma,
                nu = nu
            )
        } else {
            C <- sum(matern.covariance(
                h = s11,
                kappa = kappa, nu = nu, sigma = sigma
            ))
        }
        return(as.double(C))
    }
}


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
#' cov.approx <- matern.rational.cov(h, kappa = 10, sigma = 1, nu = 0.8, order = 2)
#' 
#' plot(h, cov.true)
#' lines(h, cov.approx, col = 2)
#'  
matern.rational.cov = function(h, 
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
    n <- length(h)
    if (nu>0 && nu<0.5)
    {
        sigma_rational = 0
        for (i in 1:length(p)){
            Sigma <- matrix(0,n,1)
            for(kk in 1:n) {
                Sigma[kk] <- matern.p(h[1],h[kk], kappa = kappa,
                                      p = p[i], alpha = alpha)
            }
            sigma_rational = sigma_rational+ r[i]*sigma^2*Sigma
        }
        Sigma <- matrix(0,n,1)
        for(kk in 1:n) {
            Sigma[kk] <- matern.k(h[1],h[kk], kappa = kappa, alpha = alpha)
        }
        sigma_rational = sigma_rational + k*sigma^2*Sigma
    }
    
    else {
        
        #k part
        Sigma <- matrix(0,n,1)
        for(kk in 1:n) {
            Sigma[kk] <- matern.k(h[1],h[kk], kappa = kappa, alpha = alpha)
        }
        sigma_rational <- k*sigma^2*Sigma
        
        # p part
        for (i in 1:length(p))
        {
            Sigma <- matrix(0,n,1)
            for(kk in 1:n) {
                Sigma[kk] <- matern.p(h[1],h[kk], kappa = kappa,
                                      p = p[i], alpha = alpha)
            }
            sigma_rational = sigma_rational+ r[i]*sigma^2*Sigma
        }
    }
    
    return(sigma_rational)
}


#Joint covariance of process and derivative for shifted Matern
#' @noRd
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

#Derivative of shifted Matern covariance 
#' @noRd
matern.k.deriv <- function(s,t,kappa,alpha, deriv = 0){
    h <- s-t
    ca <- gamma(alpha)/gamma(alpha-0.5)
    if(alpha<1){
        if(deriv == 0) {
            return((h==0)*sqrt(4*pi)*ca/kappa)    
        } else {
            stop("Not differentiable")
        }
    } else {
        fa <- floor(alpha)
        cfa <- gamma(fa)/gamma(fa-0.5)
        return(matern.derivative(h, kappa = kappa, nu = fa - 1/2, 
                                 sigma = sqrt(ca/cfa), deriv = deriv))
    }
}

#Shifted Matern covariance 
#' @noRd
matern.k <- function(s,t,kappa,alpha){
    h <- s-t
    ca <- gamma(alpha)/gamma(alpha-0.5)
    if(alpha<1){
        return((h==0)*sqrt(4*pi)*ca/kappa)
    } else {
        fa <- floor(alpha)
        cfa <- gamma(fa)/gamma(fa-0.5)
        return(matern.covariance(h, kappa = kappa, nu = fa - 1/2, 
                                 sigma = sqrt(ca/cfa)))
    }
}
