## util.R
##
##   Copyright (C) 2018, 2019, David Bolin
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Internal function to get the roots of the polynomials used in the rational approximation.
get.roots <- function(m, beta)
  {
  if (beta > 2) {
    beta <- beta - floor(beta - 1)
  }

  rb <- rep(0, m + 1)
  rc <- rep(0, m)
  if (m == 1) {
    rc    <- approx(m1table$beta, m1table$rc,beta)$y
    rb[1] <- approx(m1table$beta, m1table$rb.1,beta)$y
    rb[2] <- approx(m1table$beta, m1table$rb.2,beta)$y
    factor <- approx(m1table$beta,m1table$factor,beta)$y
  } else if (m == 2) {
    rc[1] <- approx(m2table$beta, m2table$rc.1,beta)$y
    rc[2] <- approx(m2table$beta, m2table$rc.2,beta)$y
    rb[1] <- approx(m2table$beta, m2table$rb.1,beta)$y
    rb[2] <- approx(m2table$beta, m2table$rb.2,beta)$y
    rb[3] <- approx(m2table$beta, m2table$rb.3,beta)$y
    factor <- approx(m2table$beta, m2table$factor,beta)$y
  }else if (m == 3) {
    rc[1] <- approx(m3table$beta, m3table$rc.1,beta)$y
    rc[2] <- approx(m3table$beta, m3table$rc.2,beta)$y
    rc[3] <- approx(m3table$beta, m3table$rc.3,beta)$y
    rb[1] <- approx(m3table$beta, m3table$rb.1,beta)$y
    rb[2] <- approx(m3table$beta, m3table$rb.2,beta)$y
    rb[3] <- approx(m3table$beta, m3table$rb.3,beta)$y
    rb[4] <- approx(m3table$beta, m3table$rb.4,beta)$y
    factor <- approx(m3table$beta, m3table$factor,beta)$y
  } else if (m == 4) {
    rc[1] <- approx(m4table$beta, m4table$rc.1,beta)$y
    rc[2] <- approx(m4table$beta, m4table$rc.2,beta)$y
    rc[3] <- approx(m4table$beta, m4table$rc.3,beta)$y
    rc[4] <- approx(m4table$beta, m4table$rc.4,beta)$y
    rb[1] <- approx(m4table$beta, m4table$rb.1,beta)$y
    rb[2] <- approx(m4table$beta, m4table$rb.2,beta)$y
    rb[3] <- approx(m4table$beta, m4table$rb.3,beta)$y
    rb[4] <- approx(m4table$beta, m4table$rb.4,beta)$y
    rb[5] <- approx(m4table$beta, m4table$rb.5,beta)$y
    factor <- approx(m4table$beta, m4table$factor,beta)$y
  }else{
    stop("m must be one of the values 1,2,3,4.")
  }

  return(list(rb = rb, rc = rc, factor = factor))

}

#' The Matern covariance function
#'
#' \code{matern.covariance} evaluates the Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h).}
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
#' x = seq(from = 0, to = 1, length.out = 101)
#' plot(x, matern.covariance(abs(x - 0.5), kappa = 10, nu = 1/5, sigma = 1),
#'      type = "l", ylab = "C(h)", xlab = "h")

matern.covariance <- function(h, 
                              kappa, 
                              nu, 
                              sigma)
{
  if (nu == 1/2) {
    C = sigma^2*exp(-kappa*abs(h))
  } else {
    C = (sigma^2 / (2^(nu - 1) * gamma(nu))) * ((kappa*abs(h))^nu) * besselK(kappa*abs(h), nu)
  }
  C[h == 0] = sigma^2
  return(as.matrix(C))
}

#' The 1d folded Matern covariance function
#'
#' @description 
#' \code{folded.matern.covariance.1d} evaluates the 1d folded Matern covariance function
#' over an interval \eqn{[0,L]}.
#' 
#' @details 
#' \code{folded.matern.covariance.1d} evaluates the 1d folded Matern covariance function
#' over an interval \eqn{[0,L]} under different boundary conditions. For periodic boundary conditions
#' \deqn{C_{\mathcal{P}}(h,m) = \sum_{k=-\infty}^{\infty} (C(h-m+2kL),}
#' for Neumann boundary conditions
#' \deqn{C_{\mathcal{N}}(h,m) = \sum_{k=-\infty}^{\infty} (C(h-m+2kL)+C(h+m+2kL)),}
#' and for Dirichlet boundary conditions:
#' \deqn{C_{\mathcal{D}}(h,m) = \sum_{k=-\infty}^{\infty} (C(h-m+2kL)-C(h+m+2kL)),}
#' where \eqn{C(\cdot)} is the Matern covariance function:
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h).}
#' 
#' We consider the truncation:
#' \deqn{C_{{\mathcal{P}},N}(h,m) = \sum_{k=-N}^{N} C(h-m+2kL), C_{\mathcal{N},N}(h,m) = \sum_{k=-\infty}^{\infty} (C(h-m+2kL)+C(h+m+2kL)),}
#' and
#' \deqn{C_{\mathcal{D},N}(h,m) = \sum_{k=-N}^{N} (C(h-m+2kL)-C(h+m+2kL)).}
#'
#' @param h,m Vectors of arguments of the covariance function.
#' @param kappa Range parameter.
#' @param nu Shape parameter.
#' @param sigma Standard deviation.
#' @param L The upper bound of the interval [0,L]. By default, \code{L=1}.
#' @param N The truncation parameter.
#' @param boundary The boundary condition. The possible conditions are \code{"neumann"} (default), \code{"dirichlet"} or \code{"periodic"}.
#'
#' @return A matrix with the corresponding covariance values. 
#' @export
#'
#' @examples
#' x = seq(from = 0, to = 1, length.out = 101)
#' plot(x, folded.matern.covariance.1d(rep(0.5,length(x)),x, kappa = 10, nu = 1/5, sigma = 1),
#'      type = "l", ylab = "C(h)", xlab = "h")

folded.matern.covariance.1d <- function(h, m, kappa, nu, sigma,
                                     L=1, N=10, boundary = c("neumann","dirichlet","periodic"))
{
  boundary <- tolower(boundary[1])
  if(!(boundary%in%c("neumann","dirichlet","periodic"))){
    stop("The possible boundary conditions are 'neumann', 'dirichlet' or 'periodic'!")
  }
  if(length(h)!=length(m)){
    stop("h and m should have the same length!")
  }

  s1 <- sapply(-N:N, function(j){
    h-m+2*j*L
  })
  s2 <- sapply(-N:N, function(j){
    h+m+2*j*L
  })
  if(boundary=="neumann"){
    C <- rowSums(matern.covariance(h=s1, kappa=kappa, nu=nu, sigma=sigma)+
                 matern.covariance(h=s2, kappa=kappa, nu=nu, sigma=sigma))
  } else if(boundary=="dirichlet"){
    C <- rowSums(matern.covariance(h=s1, kappa=kappa, nu=nu, sigma=sigma)-
                 matern.covariance(h=s2, kappa=kappa, nu=nu, sigma=sigma))
  } else{
    C <- rowSums(matern.covariance(h=s1, kappa=kappa, nu=nu, sigma=sigma))
  }
  
  if(length(h)==1){
    return(sum(C))
  }
  return(as.matrix(C))
}

#' The 2d folded Matern covariance function
#'
#' @description 
#' \code{folded.matern.covariance.2d} evaluates the 2d folded Matern covariance function
#' over an interval \eqn{[0,L]\times [0,L]}.
#' 
#' @details 
#' \code{folded.matern.covariance.2d} evaluates the 1d folded Matern covariance function
#' over an interval \eqn{[0,L]\times [0,L]} under different boundary conditions. For periodic boundary conditions
#' \deqn{C_{\mathcal{P}}((h_1,h_2),(m_1,m_2)) = \sum_{k_2=-\infty}^\infty \sum_{k_1=-\infty}^{\infty} (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|),}
#' for Neumann boundary conditions
#' \deqn{C_{\mathcal{N}}((h_1,h_2),(m_1,m_2)) = \sum_{k_2=-\infty}^\infty \sum_{k_1=-\infty}^{\infty} (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|)+C(\|(h_1-m_1+2k_1L,h_2+m_2+2k_2L)\|)+C(\|(h_1+m_1+2k_1L,h_2-m_2+2k_2L)\|)+C(\|(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\|)),}
#' and for Dirichlet boundary conditions:
#' \deqn{C_{\mathcal{D}}((h_1,h_2),(m_1,m_2)) = \sum_{k_2=-\infty}^\infty \sum_{k_1=-\infty}^{\infty} (C(\|(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\|)-C(\|(h_1-m_1+2k_1L,h_2+m_2+2k_2L)\|)-C(\|(h_1+m_1+2k_1L,h_2-m_2+2k_2L)\|)+C(\|(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\|)),}
#' where \eqn{C(\cdot)} is the Matern covariance function:
#' \deqn{C(h) = \frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h).}
#' 
#' We consider the truncation for \eqn{k_1,k_2} from \eqn{-N} to \eqn{N}.
#'
#' @param h,m Vectors with two coordinates.
#' @param kappa Range parameter.
#' @param nu Shape parameter.
#' @param sigma Standard deviation.
#' @param L The upper bound of the square \eqn{[0,L]\times [0,L]}. By default, \code{L=1}.
#' @param N The truncation parameter.
#' @param boundary The boundary condition. The possible conditions are \code{"neumann"} (default), \code{"dirichlet"}, \code{"periodic"} or \code{"R2"}.
#'
#' @return The correspoding covariance. 
#' @export
#'
#' @examples
#' h = c(0.5,0.5)
#' m = c(0.5,0.5)
#' folded.matern.covariance.2d(h,m, kappa = 10, nu = 1/5, sigma = 1)
#' 

folded.matern.covariance.2d <- function(h, m, kappa, nu, sigma,
                                        L=1, N=10, boundary = c("neumann","dirichlet","periodic","R2"))
{
  boundary <- tolower(boundary[1])
  if(!(boundary%in%c("neumann","dirichlet","periodic","r2"))){
    stop("The possible boundary conditions are 'neumann', 'dirichlet', 'periodic' or 'R2'!")
  }
  
  if(is.vector(h)){
    if(!is.vector(m)){
      stop("since 'h' is a vector, 'm' should be a vector!")
    }
    
    if((length(h)!=2) || (length(m)!= 2)){
      stop("The vectors h and m should have length 2!")
    }
    
  } else if(is.matrix(h)){
    if(!is.matrix(m)){
      stop("since 'h' is a matrix, 'm' should be a matrix!")
    }
    if(ncol(h)!= 2){
      stop("h must have two columns!")
    }
    if(!all(dim(h) == dim(m))){
      stop("h and m must have the same dimensions!")
    }
  }
  
  if(!is.vector(h) && !is.matrix(h)){
    stop("h should be either a vector or a matrix!")
  }
  
  list.comb <- expand.grid(-N:N, -N:N)
  
  if(is.matrix(h)){
    h_matrix_1 <- matrix(rep(h[,1],length(list.comb[,1])),nrow=nrow(h))
    h_matrix_2 <- matrix(rep(h[,2],length(list.comb[,1])),nrow=nrow(h))
    m_matrix_1 <- matrix(rep(m[,1],length(list.comb[,1])),nrow=nrow(m))
    m_matrix_2 <- matrix(rep(m[,2],length(list.comb[,1])),nrow=nrow(m))
    list_comb_matrix_1 <- t(matrix(rep(list.comb[,1], nrow(h)),ncol=nrow(h)))
    list_comb_matrix_2 <- t(matrix(rep(list.comb[,2], nrow(h)),ncol=nrow(h)))
    
    s11 <- sqrt( (h_matrix_1-m_matrix_1+2*list_comb_matrix_1*L)^2 + (h_matrix_2-m_matrix_2+2*list_comb_matrix_2*L)^2 )
    s12 <- sqrt( (h_matrix_1-m_matrix_1+2*list_comb_matrix_1*L)^2 + (h_matrix_2+m_matrix_2+2*list_comb_matrix_2*L)^2 )
    s21 <- sqrt( (h_matrix_1+m_matrix_1+2*list_comb_matrix_1*L)^2 + (h_matrix_2-m_matrix_2+2*list_comb_matrix_2*L)^2 )
    s22 <- sqrt( (h_matrix_1+m_matrix_1+2*list_comb_matrix_1*L)^2 + (h_matrix_2+m_matrix_2+2*list_comb_matrix_2*L)^2 )
    
    if(boundary=="neumann"){
      C <- rowSums(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s12, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s21, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s22, kappa=kappa, nu=nu, sigma=sigma))
    } else if(boundary=="dirichlet"){
      C <- rowSums(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma)-
                     matern.covariance(h=s12, kappa=kappa, nu=nu, sigma=sigma)-
                     matern.covariance(h=s21, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s22, kappa=kappa, nu=nu, sigma=sigma))
    } else if(boundary=="r2"){
      C <- matern.covariance(h=sqrt((h[1]-m[1])^2+(h[2]-m[2])^2), kappa = kappa, sigma = sigma, nu = nu)
    } else{
      C <- rowSums(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma))
    }
    return(as.double(C))
  } else{
    s11 <- sqrt( (h[1]-m[1]+2*list.comb[,1]*L)^2 + (h[2]-m[2]+2*list.comb[,2]*L)^2 )
    s12 <- sqrt( (h[1]-m[1]+2*list.comb[,1]*L)^2 + (h[2]+m[2]+2*list.comb[,2]*L)^2 )
    s21 <- sqrt( (h[1]+m[1]+2*list.comb[,1]*L)^2 + (h[2]-m[2]+2*list.comb[,2]*L)^2 )
    s22 <- sqrt( (h[1]+m[1]+2*list.comb[,1]*L)^2 + (h[2]+m[2]+2*list.comb[,2]*L)^2 )
    
    if(boundary=="neumann"){
      C <- sum(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s12, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s21, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s22, kappa=kappa, nu=nu, sigma=sigma))
    } else if(boundary=="dirichlet"){
      C <- sum(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma)-
                     matern.covariance(h=s12, kappa=kappa, nu=nu, sigma=sigma)-
                     matern.covariance(h=s21, kappa=kappa, nu=nu, sigma=sigma)+
                     matern.covariance(h=s22, kappa=kappa, nu=nu, sigma=sigma))
    } else if(boundary=="r2"){
      C <- matern.covariance(h=sqrt((h[1]-m[1])^2+(h[2]-m[2])^2), kappa = kappa, sigma = sigma, nu = nu)
    } else{
      C <- sum(matern.covariance(h=s11, kappa=kappa, nu=nu, sigma=sigma))
    }
    return(as.double(C))
  }

  
 

}




#' Summarise rSPDE objects
#'
#' Summary method for class "rSPDEobj"
#'
#' @param object an object of class "rSPDEobj", usually, a result of a call
#'   to \code{\link{fractional.operators}}, \code{\link{matern.operators}}, or
#'   \code{\link{spde.matern.operators}}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary rSPDEobj
summary.rSPDEobj <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.rSPDEobj"
  out$type = object$type
  if (out$type == "Matern approximation") {
    out$kappa = object$kappa
    out$sigma = object$sigma
    out$nu = object$nu
  }
  out$m = object$m
  out$n = dim(object$L)[1]
  return(out)
}

#' @param x an object of class "summary.rSPDEobj", usually, a result of a call
#'   to \code{\link{summary.rSPDEobj}}.
#' @export
#' @method print summary.rSPDEobj
#' @rdname summary.rSPDEobj
print.summary.rSPDEobj <- function(x, ...)
{

  cat("Type of approximation: ", x$type,"\n")
  if (x$type == "Matern approximation") {
    cat("Parametres of covariance function: kappa = ",
        x$kappa,", sigma = ", x$sigma, ", nu = ",x$nu, "\n")
  }
  cat("Order or rational approximation: ", x$m, "\n")
  cat("Size of discrete operators: ", x$n," x ", x$n, "\n")
}

#' @export
#' @method print rSPDEobj
#' @rdname summary.rSPDEobj
print.rSPDEobj <- function(x, ...) {
  print.summary.rSPDEobj(summary(x))
}

#' Observation matrix for finite element discretization on R
#'
#' A finite element discretization on R can be written as
#' \eqn{u(s) = \sum_i^n u_i \varphi_i(s)}{u(s) = \sum_i^n u_i \varphi_i(s)}
#' where \eqn{\varphi_i(s)} is a piecewise linear "hat function" centered at location
#' \eqn{x_i}{x_i}. This function computes an \eqn{m\times n}{m x n} matrix \eqn{A}{A}
#' that links the basis function in the expansion to specified locations
#' \eqn{s = (s_1,\ldots, s_m)} in the domain through
#' \eqn{A_ij = \varphi_j(s_i)}{A_ij = \varphi_j(s_i)}.
#'
#' @param x The locations of the nodes in the FEM discretization.
#' @param loc The locations \eqn{(s_1,\ldots, s_m)}
#'
#' @return The sparse matrix \code{A}.
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso \code{\link{rSPDE.fem1d}}
#'
#' @examples
#' #create mass and stiffness matrices for a FEM discretization on [0,1]
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #create the observation matrix for some locations in the domain
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
rSPDE.A1d <- function(x, loc)
{
  if (min(loc) < min(x) || max(loc) > max(x))
    stop("locations outside support of basis")

  n.x  <- length(x)
  n.loc <- length(loc)
  i <- as.vector(cbind(1:n.loc, 1:n.loc))
  j <- matrix(0, n.loc, 2)
  vals <- matrix(1, n.loc, 2)
  for (ii in seq_len(n.loc)) {
    j[ii,1] <- sum(sum((loc[ii] - x) >= 0))
    vals[ii,1] <- loc[ii] - x[j[ii,1]]
    j[ii,2] <- j[ii,1] + 1
    if (j[ii,2] <= n.x) {
      vals[ii,2] <- x[j[ii,2]] - loc[ii]
    } else {
      j[ii,2] = j[ii,2] - 2
    }
  }
  j <- as.vector(j)
  vals <- as.vector(matrix(1 - vals/rowSums(vals)))

  A <- sparseMatrix(i = i, j = j, x = vals, dims = c(n.loc, n.x))
  return(A)
}


#' Finite element calculations for problems on R
#'
#' This function computes mass and stiffness matrices for a FEM approximation on R, assuming 
#' Neumann boundary conditions.
#' These matrices are needed when discretizing the operators in rational approximations.
#'
#' @param x Locations of the nodes in the FEM approximation.
#'
#' @return The function returns a list with the following elements
#' \item{G }{The stiffness matrix.}
#' \item{C }{The mass matrix.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso \code{\link{rSPDE.A1d}}
#' @examples
#' #create mass and stiffness matrices for a FEM discretization on [0,1]
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
rSPDE.fem1d <- function(x)
{
  n = length(x)
  d <- c(Inf, diff(x))
  dm1 = c(d[2:n], Inf)
  G = -bandSparse(n = n, m = n, k = c(-1, 0, 1),
                  diagonals = cbind(1/dm1, -(1/dm1 + 1/d), 1/dm1))
  C = bandSparse(n = n, m = n, k = c(-1, 0, 1),
                 diagonals = cbind(dm1/6, (dm1 + d)/3, d/6))
  C[1, 1:2] <- c(d[2], d[2]/2)/3
  C[n, (n - 1):n] <- c(d[n]/2, d[n])/3

  return(list(G = G, C = C))
}

#' Warnings free loading of add-on packages
#'
#' Turn off all warnings for require(), to allow clean completion
#' of examples that require unavailable Suggested packages.
#'
#' @param package The name of a package, given as a character string.
#' @param lib.loc a character vector describing the location of R library trees
#' to search through, or \code{NULL}.  The default value of \code{NULL}
#' corresponds to all libraries currently known to \code{.libPaths()}.
#' Non-existent library trees are silently ignored.
#' @param character.only a logical indicating whether \code{package} can be
#' assumed to be a character string.
#'
#' @return \code{require.nowarnings} returns (invisibly) \code{TRUE} if it succeeds, otherwise \code{FALSE}
#' @details \code{require(package)} acts the same as
#' \code{require(package, quietly = TRUE)} but with warnings turned off.
#' In particular, no warning or error is given if the package is unavailable.
#' Most cases should use \code{requireNamespace(package, quietly = TRUE)} instead,
#' which doesn't produce warnings.
#' @seealso \code{\link{require}}
#' @export
#' @examples
#' ## This should produce no output:
#' if (require.nowarnings(nonexistent)) {
#'   message("Package loaded successfully")
#' }

require.nowarnings <- function(package, lib.loc = NULL, character.only = FALSE)
{
  if (!character.only)
    package <- as.character(substitute(package))
  suppressWarnings(
    require(package,
            lib.loc = lib.loc,
            quietly = TRUE,
            character.only = TRUE)
  )
}

#' @name get.inital.values.rSPDE
#' @title Initial values for log-likelihood optimization in rSPDE models
#' with a latent stationary Gaussian Matern model
#' @description Auxiliar function to obtain domain-based initial values for
#' log-likelihood optimization in rSPDE models
#' with a latent stationary Gaussian Matern model
#' @param mesh An in INLA mesh
#' @param mesh.range The range of the mesh. 
#' @param dim The dimension of the domain.
#' @param include.nu Should we also provide an initial guess for nu?
#' @param log.scale Should the results be provided in log scale?
#' @param nu_upper_bound Should an upper bound for nu be considered?
#' @param include.tau Should tau be returned instead of sigma?
#' @return A vector of the form (theta_1,theta_2,theta_3) or where
#' theta_1 is the initial guess for tau, theta_2 is the initial guess for kappa
#' and theta_3 is the initial guess for nu.
#' @export
#' 

get.inital.values.rSPDE <- function(mesh = NULL, mesh.range = NULL, 
                                    dim = NULL,
                                    include.nu = TRUE, log.scale = TRUE,
                                    include.tau = FALSE,
                                    nu_upper_bound = NULL){
  if(is.null(mesh)&&is.null(mesh.range)){
    stop("You should either provide mesh or mesh.range!")
  }
  
  if(is.null(mesh)&&is.null(dim)){
    stop("If you don't provide mesh, you have to provide dim!")
  }
  if(include.nu){
    if(!is.null(nu_upper_bound)){
      nu <- min(1,nu_upper_bound/2)
    } else{
      nu <- 1
    }
  }
  
  if(!is.null(mesh)){
    dim <- get_inla_mesh_dimension(inla_mesh = mesh)
    mesh.range = ifelse(dim == 2, (max(c(diff(range(mesh$loc[,
                                                           1])), diff(range(mesh$loc[, 2])), diff(range(mesh$loc[,
                                                                                                                 3]))))), diff(mesh$interval))
  }
  
  range.nominal <- mesh.range * 0.2
    
  kappa <- sqrt(8 * nu)/range.nominal
  
  if(include.tau){
    tau <- sqrt(gamma(nu) / (kappa^(2*nu) * (4*pi)^(dim /2) * gamma(nu + dim/2)))
    
    initial <- c(tau, kappa, nu)
  } else{
    initial <- c(1, kappa, nu)
  }  

  if(log.scale){
    return(log(initial))
  } else{
    return(initial)
  }
}




#' @name cut_decimals
#' @title Approximation function for covariance-based rSPDE models
#' @description Approximation function to be used to compute the precision matrix for covariance-based rSPDE models
#' @param nu A real number
#' @return An approximation
#' @noRd

cut_decimals <- function(nu) {
  temp <- nu - floor(nu)
  if (temp < 10 ^ (-3)) {
    temp <- 10 ^ (-3)
  }
  if (temp > 0.999) {
    temp <- 0.999
  }
  return(temp)
}

#' @name check_class_inla_rspde
#' @title Check if the object inherits from inla.rspde class
#' @description Check if the object inherits from inla.rspde class
#' @param model A model to test if it inherits from inla.rspde
#' @return Gives an error if the object does not inherit from inla.rspde
#' @noRd

check_class_inla_rspde <- function(model){
  if(!inherits(model, "inla.rspde")){
    stop("You should provide a rSPDE model!")
  }
}

#' @name get_inla_mesh_dimension
#' @title Get the dimension of an INLA mesh
#' @description Get the dimension of an INLA mesh
#' @param inla_mesh An INLA mesh
#' @return The dimension of an INLA mesh.
#' @noRd


get_inla_mesh_dimension <- function(inla_mesh){
  cond1 <- inherits(inla_mesh,"inla.mesh.1d")
  cond2 <- inherits(inla_mesh,"inla.mesh")
  stopifnot(cond1 || cond2)
  if(inla_mesh$manifold == "R1"){
      d = 1
  } else if(inla_mesh$manifold == "R2"){
      d = 2
  } else{
    stop("The mesh should be from a flat manifold.")
  }
  return(d)
}


#' @name fem_mesh_order_1d
#' @title Get fem_mesh_matrices for 1d inla.mesh objects
#' @description Get fem_mesh_matrices for 1d inla.mesh objects
#' @param inla_mesh An INLA mesh
#' @param m_order the order of the FEM matrices
#' @return A list with fem_mesh_matrices
#' @noRd


fem_mesh_order_1d <- function(inla_mesh, m_order){
  fem_mesh <- rSPDE.fem1d(inla_mesh[["loc"]])
  C <- fem_mesh$C
  C <- Matrix::Diagonal(dim(C)[1], rowSums(C))
  C <- INLA::inla.as.sparse(C)
  G <- fem_mesh$G
  Gk <- list()
  Ci <- C
  Ci@x <- 1/(C@x)
  
  GCi <- G%*%Ci
  Gk[[1]] <- G
  # determine how many G_k matrices we want to create
  if(m_order>1){
    for (i in 2:m_order){
      Gk[[i]] <- GCi %*% Gk[[i-1]]
    }
  }
  
  # create a list contains all the finite element related matrices
  fem_mesh_matrices <- list()
  fem_mesh_matrices[["c0"]] <- C
  
  for(i in 1:m_order){
    fem_mesh_matrices[[paste0("g",i)]] <- Gk[[i]]
  }
  return(fem_mesh_matrices)
}

#' @name get.sparsity.graph.rspde
#' @title Sparsity graph for rSPDE models
#' @description Creates the sparsity graph for rSPDE models
#' @param mesh An INLA mesh, optional
#' @param fem_mesh_matrices A list containing the FEM-related matrices. The list should contain elements C, G, G_2, G_3, etc. Optional,
#' should be provided if mesh is not provided.
#' @param dim The dimension, optional. Should be provided if mesh is not provided.
#' @param nu The smoothness parameter
#' @param force_non_integer Should nu be treated as non_integer?
#' @param rspde_order The order of the covariance-based rational SPDE approach.
#' @param sharp The graph should have the correct sparsity (costs more to perform
#' a sparsity analysis) or an upper bound for the sparsity?
#' @return The sparsity graph for rSPDE models to be used in R-INLA interface.
#' @export

get.sparsity.graph.rspde <- function(mesh=NULL,
                                     fem_mesh_matrices=NULL,
                                     nu,
                                     force_non_integer = FALSE,
                                     rspde_order = 2,
                                     sharp = TRUE,
                                     dim=NULL){
  
  if(!is.null(mesh)){
    stopifnot(inherits(mesh,"inla.mesh"))
    if(mesh$manifold == "R1"){
      dim = 1
    } else if(mesh$manifold == "R2"){
      dim = 2
    } else{
      stop("The mesh should be from a flat manifold.")
    }
  } else if(is.null(dim)){
    stop("If an INLA mesh is not provided, you should provide the dimension!")
  }
  
  alpha = nu + dim/2
  
  m_alpha = floor(alpha)
  
  integer_alpha <- (alpha%%1 == 0)
  
  if(force_non_integer){
    integer_alpha = FALSE
  }
  
  if(!is.null(fem_mesh_matrices)){
    if(integer_alpha){
      return(fem_mesh_matrices[[paste0("g",m_alpha)]])
    } else{
      if(sharp){
        if(m_alpha>0){
          return(bdiag(kronecker(diag(rep(1,rspde_order)),
                                 fem_mesh_matrices[[paste0("g",m_alpha+1)]]), 
                       fem_mesh_matrices[[paste0("g",m_alpha)]]))
        } else{
          return(bdiag(kronecker(diag(rep(1,rspde_order)),
                                 fem_mesh_matrices[["g1"]]), 
                       fem_mesh_matrices[["c0"]]))
        }

      } else{
        return(kronecker(diag(rep(1,rspde_order+1)),
                         fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      }
    }
  } else if(!is.null(mesh)){
    if(integer_alpha){
      fem_mesh_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha)
      return(fem_mesh_matrices[[paste0("g",m_alpha)]])
    } else{
      if(dim == 2){
        fem_mesh_matrices <- INLA::inla.mesh.fem(mesh, order = m_alpha+1)
      } else{
        fem_mesh_matrices <- fem_mesh_order_1d(mesh, m_order = m_alpha+1)
      }

      
      if(sharp){
        if(m_alpha>0){
          return(bdiag(kronecker(diag(rep(1,rspde_order)),
                                 fem_mesh_matrices[[paste0("g",m_alpha+1)]]), 
                       fem_mesh_matrices[[paste0("g",m_alpha)]]))
        } else{
          return(bdiag(kronecker(diag(rep(1,rspde_order)),
                                 fem_mesh_matrices[["g1"]]), 
                       fem_mesh_matrices[["c0"]]))
        }

      } else{
        return(kronecker(diag(rep(1,rspde_order+1)),
                         fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      }
    } 
  } else {
    stop("You should provide either mesh or fem_mesh_matrices!")
  } 
}


#' @name build_sparse_matrix_rspde
#' @title Create sparse matrix from entries and graph
#' @description Create sparse matrix from entries and graph
#' @param entries The entries of the precision matrix
#' @param graph The sparsity graph of the precision matrix
#' @return index for rSPDE models.
#' @noRd

build_sparse_matrix_rspde <- function(entries, graph){
  if(!is.null(graph)){
    graph = as(graph, "dgTMatrix")
    idx <- which(graph@i <= graph@j)
    Q = Matrix::sparseMatrix(i=graph@i[idx], j = graph@j[idx], x= entries,
                             symmetric=TRUE, index1 = FALSE)
  }
  return(Q)
}


#' @name analyze_sparsity_rspde
#' @title Analyze sparsity of matrices in the rSPDE approach
#' @description Auxiliar function to analyze sparsity of matrices in the rSPDE approach
#' @param nu_upper_bound Upper bound for the smoothness parameter
#' @param dim The dimension of the domain
#' @param rspde_order The order of the rational approximation
#' @param fem_mesh_matrices A list containing FEM-related matrices. The list should contain elements c0, g1, g2, g3, etc.
#' @param include_lower_order Logical. Should the lower-order terms be included? They are needed for the cases
#' when alpha = nu + d/2 is integer or for when sharp is set to TRUE.
#' @param include_higher_order Logical. Should be included for when nu is estimated or for when alpha = nu + d/2 is not an integer.
#' @return A list containing informations on sparsity of the precision matrices
#' @noRd

analyze_sparsity_rspde <- function(nu_upper_bound, dim, rspde_order,
                                   fem_mesh_matrices,
                                   include_lower_order = TRUE, 
                                   include_higher_order = TRUE){
  beta = nu_upper_bound / 2 + dim / 4
  
  m_alpha = floor(2 * beta)
  
  positions_matrices <- list()
  
  C_list <- symmetric_part_matrix(fem_mesh_matrices$c0)
  G_1_list <- symmetric_part_matrix(fem_mesh_matrices$g1)
  if(m_alpha<2){
      G_2_list <- symmetric_part_matrix(fem_mesh_matrices[["g2"]])
    }
  if(m_alpha>1){
    for(j in 2:(m_alpha)){
      assign(paste0("G_",j,"_list"), symmetric_part_matrix(fem_mesh_matrices[[paste0("g",j)]]))
    }
    
    if(include_higher_order){
      assign(paste0("G_",m_alpha+1,"_list"), symmetric_part_matrix(fem_mesh_matrices[[paste0("g",m_alpha+1)]]))
      
      positions_matrices[[1]] <-  match(C_list$M, get(paste0("G_",m_alpha+1,"_list"))[["M"]])
    }
  }
  
  idx_matrices <- list()
  
  idx_matrices[[1]] = C_list$idx
  
  if(m_alpha>0){
  for(i in 1:m_alpha){
    if(include_higher_order){
      positions_matrices[[i+1]] <- match(get(paste0("G_",i,"_list"))[["M"]], get(paste0("G_",m_alpha+1,"_list"))[["M"]])
    }
    idx_matrices[[i+1]] = get(paste0("G_",i,"_list"))[["idx"]]
  }
  }
  
  if(include_higher_order){
    idx_matrices[[m_alpha+2]] = get(paste0("G_",m_alpha+1,"_list"))[["idx"]]
  }
  
  if(include_lower_order){
    positions_matrices_less <- list()
    if(m_alpha>0){
      positions_matrices_less[[1]] <-  match(C_list$M, get(paste0("G_",m_alpha,"_list"))[["M"]])
    } else{
      positions_matrices_less[[1]] <-  match(C_list$M, get(paste0("G_",1,"_list"))[["M"]])
    }

    if(m_alpha > 1){
      for(i in 1:(m_alpha-1)){
        positions_matrices_less[[i+1]] <- match(get(paste0("G_",i,"_list"))[["M"]], get(paste0("G_",m_alpha,"_list"))[["M"]])
      }
    } else if (m_alpha == 1){
      positions_matrices_less[[2]] <-  1:length(get(paste0("G_",m_alpha,"_list"))[["M"]])
    }
  } else{
    positions_matrices_less <- NULL
  }
  
  return(list(positions_matrices = positions_matrices,
              idx_matrices = idx_matrices,
              positions_matrices_less = positions_matrices_less))
}

#' @name symmetric_part_matrix
#' @title Gets the upper triangular part of a matrix
#' @description Gets the upper triangular part of a matrix
#' @param M A matrix or a sparse matrix
#' @return A sparse matrix formed by the upper triangular part of \code{M}.
#' @noRd

symmetric_part_matrix <- function(M){
  M  <- as(M, "dgTMatrix")
  idx <- which(M@i <= M@j)
  sM <- cbind(M@i[idx], M@j[idx])
  colnames(sM) <- NULL
  return(list(M = split(sM, seq(nrow(sM))), idx = idx))
}


#' @name create_summary_from_density
#' @title Creates a summary from a density data frame
#' @description Auxiliar function to create summaries from density data drames
#' @param density_df A density data frame
#' @param name Name of the parameter
#' @return A data frame containing a basic summary
#' @noRd

create_summary_from_density <- function(density_df, name){
  min_x <- min(density_df[,"x"])
  max_x <- max(density_df[,"x"])
  denstemp <- function(x){
    dens <- sapply(x, function(z){
      if(z<min_x){
        return(0)
      } else if (z>max_x){
        return(0)
      } else{
        return(approx(x = density_df[,"x"], y = density_df[,"y"], xout = z)$y)
      }
    })
    return(dens)
  } 
  
  ptemp <- function(q){
    prob_temp <- sapply(q, function(v){
      if(v<=min_x){
        return(0)
      } else if(v>=max_x){
        return(1)
      } else{
        stats::integrate(f = denstemp,lower = min_x, upper = v,
                         subdivisions = min(nrow(density_df),500))$value
      }
    })
    return(prob_temp)
  }
  
  mean_temp <- stats::integrate(f = function(z){denstemp(z)*z},lower = min_x, upper = max_x,
                                subdivisions = nrow(density_df))$value
  
  sd_temp <- sqrt(stats::integrate(f = function(z){denstemp(z)*(z-mean_temp)^2},lower = min_x, upper = max_x,
                                   subdivisions = nrow(density_df))$value)
  
  mode_temp <- density_df[which.max(density_df[,"y"]),"x"]
  
  qtemp <- function(p){
    quant_temp <- sapply(p, function(x){
      if(x < 0 | x > 1){
        warn_qtemp <- TRUE
        return(NaN)} else{
          return(stats::uniroot(function(y) {
            ptemp(y) - x
          },lower = min_x, upper = max_x)$root)
        }
    }
    )
    return(quant_temp)
  }
  
  out <- data.frame(mean = mean_temp, sd = sd_temp, `0.025quant`=qtemp(0.025), 
                    `0.5quant`=qtemp(0.5), `0.975quant`=qtemp(0.975), mode = mode_temp)
  rownames(out) <- name
  colnames(out) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
  return(out)
}






#' Summarise CBrSPDE objects
#'
#' Summary method for class "CBrSPDEobj"
#'
#' @param object an object of class "CBrSPDEobj", usually, a result of a call
#'   to \code{\link{matern.operators}}.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary CBrSPDEobj
#' @examples
#' #Compute the covariance-based rational approximation of a 
#' #Gaussian process with a Matern covariance function on R
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation of covariance function at 0.5
#' tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(1/2) * gamma(nu+1/2)))
#' op_cov <- matern.operators(C=fem$C, G=fem$G,nu=nu,
#' kappa=kappa,sigma=sigma,d=1,m=2)
#' 
#' op_cov
summary.CBrSPDEobj <- function(object, ...)
{
  out <- list()
  class(out) <- "summary.CBrSPDEobj"
  out$type = object$type
  out$kappa = object$kappa
  out$sigma = object$sigma
  out$nu = object$nu
  out$m = object$m
  out$n = dim(object$C)[1]
  out[["type_rational_approximation"]] = object[["type_rational_approximation"]]
  return(out)
}

#' @param x an object of class "summary.CBrSPDEobj", usually, a result of a call
#'   to \code{\link{summary.CBrSPDEobj}}.
#' @export
#' @method print summary.CBrSPDEobj
#' @rdname summary.CBrSPDEobj
print.summary.CBrSPDEobj <- function(x, ...)
{
  
  cat("Type of approximation: ", x$type,"\n")
  cat("Type of rational approximation: ", x[["type_rational_approximation"]],"\n")
  cat("Parameters of covariance function: kappa = ",
        x$kappa,", sigma = ", x$sigma, ", nu = ",x$nu, "\n")
  cat("Order or rational approximation: ", x$m, "\n")
  cat("Size of discrete operators: ", x$n," x ", x$n, "\n")
}

#' @export
#' @method print CBrSPDEobj
#' @rdname summary.CBrSPDEobj
print.CBrSPDEobj <- function(x, ...) {
  print.summary.CBrSPDEobj(summary(x))
}



#' @name get_rational_coefficients
#' @title Get matrix with rational coefficients
#' @description Get matrix with rational coefficients
#' @param order order of the rational approximation
#' @param type_rational_approx Type of the rational approximation. Options are "chebfun", "brasil" and "chebfunLB"
#' @return A matrix with rational approximations.
#' @noRd

get_rational_coefficients <- function(order, type_rational_approx){
  if(type_rational_approx == "chebfun"){
    mt <- get(paste0("m", order, "t"))
  } else if(type_rational_approx == "brasil"){
    mt <- get(paste0("m_brasil", order, "t"))
  } else if(type_rational_approx == "chebfunLB"){
    mt <- get(paste0("m_chebfun", order, "t"))
  } else{
    stop("The options are 'chebfun', 'brasil' and 'chebfunLB'!")
  }
  return(mt)
}




#' Changing the type of the rational approximation
#'
#' @param x A \code{CBrSPDE} or an \code{rpsde.inla} object
#' @param value The type of rational approximation. The current options are "chebfun", "brasil" and "chebfunLB"
#'
#' @return An object of the same class with the new rational approximation.
#' @export
#'
`rational.type<-` <- function(x, value){
  object <- x
  
  type_rational_approximation <- value
  type_rational_approximation <- type_rational_approximation[[1]]
  if(!(type_rational_approximation %in% c("chebfun","brasil","chebfunLB"))){
    stop('The possible types are "chebfun", "brasil" and "chebfunLB"!')
  }
  if(inherits(object, "CBrSPDEobj")){
    model <- update(x, type_rational_approximation = value)
  } else if(inherits(object, "inla.rspde")){
    
    fem_mesh <- object$fem_mesh
    optimize <- object$optimize
    rgeneric_type <- object$rgeneric_type
    d <- object$dim
    est_nu <- object$est_nu
    n.spde <- object$n.spde
    nu_upper_bound <- object$nu_upper_bound
    prior.nu.dist <- object$prior.nu.dist
    graph_opt <- object$f$rgeneric$definition(cmd="graph")
    
    nu <- object[["nu"]]
    prior.kappa <- object$prior.kappa
    prior.nu <- object$prior.nu 
    prior.tau <- object$prior.tau 
    start.lkappa <- object$start.lkappa 
    start.ltau <- object$start.ltau
    start.nu <- object$start.nu 
    rspde_order <- object$rspde_order
    sharp <- object$sharp
    fem_matrices <- object$fem_matrices
    C <- fem_mesh$c0
    debug <- object$debug
    
    if(rgeneric_type == "general"){
      model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_general,
                                          nu_upper_bound=nu_upper_bound,
                                          fem_matrices=fem_matrices, 
                                          graph_opt=graph_opt,
                                          sharp=sharp,
                                          prior.kappa=prior.kappa,
                                          prior.nu=prior.nu,
                                          prior.tau=prior.tau,
                                          start.lkappa = start.lkappa,
                                          start.nu = start.nu,
                                          start.ltau = start.ltau,
                                          type.rational.approx=type_rational_approximation,
                                          d = d, rspde_order = rspde_order, 
                                          prior.nu.dist=prior.nu.dist,
                                          n=ncol(C)*(rspde_order+1), 
                                          debug=debug,
                                          do_optimize=optimize, optimize=optimize)
    } else if(rgeneric_type == "frac_alpha"){
      model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_frac_alpha,
                                          nu = nu,
                                          fem_matrices=fem_matrices, 
                                          graph_opt=graph_opt,
                                          sharp=sharp,
                                          prior.kappa=prior.kappa,
                                          prior.nu=prior.nu,
                                          prior.tau=prior.tau,
                                          start.lkappa = start.lkappa,
                                          start.ltau = start.ltau,
                                          type.rational.approx=type.rational.approx,
                                          d = d, rspde_order = rspde_order, 
                                          n=ncol(C)*(rspde_order+1), 
                                          debug=debug,
                                          do_optimize=optimize, optimize=optimize)
    } else if(rgeneric_type == "int_alpha"){
      model <- INLA::inla.rgeneric.define(inla.rgeneric.cov_rspde_int_alpha,
                                          nu = nu,
                                          fem_matrices=fem_matrices, 
                                          graph_opt=graph_opt,
                                          prior.kappa=prior.kappa,
                                          prior.nu=prior.nu,
                                          prior.tau=prior.tau,
                                          type.rational.approx=type.rational.approx,
                                          start.lkappa = start.lkappa,
                                          start.ltau = start.ltau,
                                          d = d,
                                          n=ncol(C), 
                                          debug=debug,
                                          do_optimize=optimize, optimize=optimize)
    } else{
      stop("Something is wrong with the rspde model!")
    }
    model$optimize = optimize
    model$nu <- nu
    model$prior.kappa <- prior.kappa
    model$prior.nu <- prior.nu
    model$prior.tau <- prior.tau
    model$start.lkappa <- start.lkappa
    model$start.ltau <- start.ltau
    model$start.nu <- start.nu
    model$rspde_order <- rspde_order
    class(model) <- c(class(model), "inla.rspde")
    model$dim = d
    model$est_nu = est_nu
    model$n.spde = n.spde
    model$nu_upper_bound = nu_upper_bound
    model$prior.nu.dist = prior.nu.dist
    model$type.rational.approx = type_rational_approximation
    model$sharp = sharp
    model$debug = debug
    model$mesh = object$mesh
    model$fem_mesh <- fem_mesh
    model$fem_matrices <- fem_matrices
  } else{
    stop("The object must be of class 'CBrSPDE' or 'inla.rspde'!")
  }
  return(model)
}



#' Get type of rational approximation.
#'
#' @param object A \code{CBrSPDEobj} object or an \code{inla.rspde} object.
#'
#' @return The type of rational approximation.
#' @export
#'
rational.type <- function(object){
    if(inherits(object, "CBrSPDEobj")){
      return(object$type_rational_approximation)
    } else if(inherits(object, "inla.rspde")){
      return(object$type.rational.approx)
    } else if(inherits(object, "rSPDEobj")){
      return("chebfun")
    } else {
      stop("Not a valid rSPDE object!")
    }
}


#' Changing the order of the rational approximation
#'
#' @param x A \code{CBrSPDE} or an \code{rpsde.inla} object
#' @param value The order of rational approximation.
#' 
#' @return An object of the same class with the new order of rational approximation.
#' @export
#'
`rational.order<-` <- function(x, value){
  object <- x

  rspde_order <- value
  rspde_order <- rspde_order[[1]]

  if(inherits(object, "CBrSPDEobj") || inherits(object, "rSPDEobj")){
    model <- update(object, user_m = rspde_order)
  } else if(inherits(object, "inla.rspde")){
    
    optimize <- object$optimize
    nu_upper_bound <- object$nu_upper_bound
    prior.nu.dist <- object$prior.nu.dist
    mesh = object$mesh
    nu <- object[["nu"]]
    prior.kappa <- object$prior.kappa
    prior.nu <- object$prior.nu 
    prior.tau <- object$prior.tau 
    start.lkappa <- object$start.lkappa 
    start.ltau <- object$start.ltau
    start.nu <- object$start.nu 
    type_rational_approximation <- object$type.rational.approx
    sharp <- object$sharp
    debug <- object$debug
    
    
    model <- rspde.matern(mesh, 
      nu_upper_bound = nu_upper_bound,
      rspde_order = rspde_order,
      nu = nu, sharp = sharp,
      debug = debug,
      optimize=optimize, 
      prior.kappa = prior.kappa,
      prior.nu = prior.nu,
      prior.tau = prior.tau,
      start.lkappa = start.lkappa,
      start.nu = start.nu,
      start.ltau = start.ltau,
      prior.nu.dist = prior.nu.dist,
      type.rational.approx = type_rational_approximation
    )
    
  } else if(!is.null(attr(object, "inla.rspde.Amatrix"))){
      n_temp <- ncol(object)
      old_rspde_order <- attr(object,"rspde_order")
      orig_dim <- n_temp/(old_rspde_order+1)
      A <- object[, 1:orig_dim]
      Abar = kronecker(matrix(1,1,rspde_order+1),A)
      attr(Abar, "inla.rspde.Amatrix") <- TRUE
      attr(Abar,"rspde_order") <- rspde_order
      integer_nu <- attr(object, "integer_nu")
      if(integer_nu && rspde_order > 0){
        warning("The order was not changed since there is no rational approximation (an integer model was considered).")
        return(object)
      }
      attr(Abar, "integer_nu") <- integer_nu
      return(Abar)
  } else if(inherits(object, "inla.rspde.index")) {

    integer_nu <- attr(object, "integer_nu") 
    
    if(integer_nu && rspde_order > 0){
      warning("The order was not changed since there is no rational approximation (an integer model was considered).")
      return(object)
    }
    
    n_mesh <- attr(object, "n.mesh")
    name <- attr(object, "name")
    n.group <- attr(object, "n.group")
    n.repl <- attr(object, "n.repl")
    
    factor_rspde = rspde_order+1
    
    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    
    out <- list()
    out[[name]] <- as.vector(sapply(1:factor_rspde, function(i){rep(rep(((i-1)*n_mesh+1):(i*n_mesh), times = n.group), times = n.repl)}))
    out[[name.group]] <- rep(rep(rep(1:n.group, each = n_mesh), times = n.repl), times = factor_rspde)
    out[[name.repl]] <- rep(rep(1:n.repl, each = n_mesh * n.group), times = factor_rspde)
    class(out) <- c(class(out), "inla.rspde.index")
    attr(out, "rspde_order") <- rspde_order
    attr(out, "integer_nu") <- integer_nu
    attr(out, "n.mesh") <- n_mesh
    attr(out, "name") <- name
    attr(out, "n.group") <- n.group
    attr(out, "n.repl") <- n.repl
    return(out)
    
  } else{
    stop("The object must be of class 'CBrSPDE' or 'inla.rspde'!")
  }
  return(model)
}


#' Get the order of rational approximation.
#'
#' @param object A \code{CBrSPDEobj} object or an \code{inla.rspde} object.
#'
#' @return The order of rational approximation.
#' @export
#'
rational.order <- function(object){
  if(inherits(object, "CBrSPDEobj") || inherits(object, "rSPDEobj")){
    return(object$m)
  } else if(inherits(object, "inla.rspde")){
    return(object$rspde_order)
  } else if(!is.null(attr(object, "inla.rspde.Amatrix"))){
    return(attr(object,"rspde_order"))
  } else if(inherits(object, "inla.rspde.index")){
    return(attr(object,"rspde_order"))
  } else{
    stop("Not a valid rSPDE object!")
  }
}
