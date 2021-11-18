## fractional.computations.R
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

#' @rdname simulate.rSPDEobj
#' @export
simulate <- function(object, nsim,...) {
  UseMethod("simulate", object)
}

#' Simulation of a fractional SPDE using a rational SPDE approximation
#'
#' The function samples a Gaussian random field based on a pre-computed rational SPDE approximation.
#'
#' @param object The rational SPDE approximation, computed using \code{\link{fractional.operators}},
#' \code{\link{matern.operators}}, or \code{\link{spde.matern.operators}}.
#' @param nsim The number of simulations.
#' @param ... Currently not used.
#'
#' @return A matrix with the \code{n} samples as columns.
#' @export
#' @method simulate rSPDEobj
#'
#' @examples
#' #Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation
#' op <- matern.operators(kappa = kappa, sigma = sigma,
#'                        nu = nu, G=fem$G, C=fem$C, d = 1)
#'
#' #Sample the model and plot the result
#' Y <- simulate(op)
#' plot(x, Y, type = "l", ylab = "u(x)", xlab = "x")

simulate.rSPDEobj <- function(object, nsim = 1 ,...)
{
  if (class(object) != "rSPDEobj")
    stop("input op is not of class rSPDEobj")
  m <- dim(object$Q)[1]
  z <- rnorm(nsim * m)
  dim(z) <- c(m, nsim)
  x <- Qsqrt.solve(object,z)
  x <- Pr.mult(object,x)
  
  return(x)
}


#' @name update.cov_rSPDEobj
#' @title Update parameters of cov_rSPDE objects
#' @description Function to change the parameters of a cov_rSPDE object
#' @param object A cov_rSPDE object
#' @param user_kappa If non-null, update the range parameter of the covariance function.
#' @param user_tau If non-null, update the standard deviation of the covariance function.
#' @param user_nu If non-null, update the shape parameter of the covariance function.
#' @param user_m If non-null, update the order of the rational approximation, which needs to be a positive integer.
#' @param ... Currently not used.
#' @return A cov_rSPDE object
#' @method update cov_rSPDEobj
#' @export


update.cov_rSPDEobj <- function(object, user_nu = NULL,
                                user_kappa = NULL,
                                user_tau = NULL,
                                user_m = NULL, ...){
  new_object <- object
  d <- object$d
  
  ## get parameters
  if(!is.null(user_nu)){
    if(!is.numeric(user_nu)){
      stop("nu should be a number!")
    }
    if(length(user_nu)>1){
      stop("nu should be a number!")      
    }
    new_object$nu <- user_nu
    nu <- user_nu
    alpha <- nu + d/2
    m_alpha <- max(1,floor(alpha))
    m_order <- ifelse(alpha%%1==0, m_alpha, m_alpha+1)
    
    fem_mesh_matrices <- object$fem_mesh_matrices
    
    if( m_order + 1 > length(object$fem_mesh_matrices) ){
      old_m_order <- length(object$fem_mesh_matrices) - 1
      GCi <- object$GCi
      for (i in (old_m_order+1):m_order){
        fem_mesh_matrices[[paste0("g",i)]] <- GCi %*% fem_mesh_matrices[[paste0("g",i-1)]]
      }
    }
  }
  
  new_object[["fem_mesh_matrices"]] <- fem_mesh_matrices 
  
  
  if(!is.null(user_kappa)){
    if(!is.numeric(user_kappa)){
      stop("kappa should be a number!")
    }
    if(length(user_kappa)>1){
      stop("kappa should be a number!")      
    }
    new_object$kappa <- user_kappa
  }
  
  if(!is.null(user_tau)){
    if(!is.numeric(user_tau)){
      stop("tau should be a number!")
    }
    if(length(user_tau)>1){
      stop("tau should be a number!")      
    }
    new_object$tau <- user_tau
  }
  
  
  if(!is.null(user_m)){
    if(!is.numeric(user_m)){
      stop("m should be a number greater or equal to 1")
    }
    if(length(user_m)>1){
      stop("m should be a number greater or equal to 1")
    }
    if(user_m < 1){
      stop("m should be a number greater or equal to 1")
    }
    new_object$m <- as.integer(user_m)
  }
  
  return(new_object)
}

#' @name simulate.cov_rSPDEobj
#' @title Simulation of a fractional SPDE using a rational SPDE approximation
#' @description The function samples a Gaussian random field based on a pre-computed rational SPDE approximation.
#' @param object A cov_rSPDE object
#' @param nsim The number of simulations.
#' @param user_kappa If non-null, update the range parameter of the covariance function.
#' @param user_tau If non-null, update the standard deviation of the covariance function.
#' @param user_nu If non-null, update the shape parameter of the covariance function.
#' @param user_m If non-null, update the order of the rational approximation, which needs to be a positive integer.
#' @param pivot Should pivoting be used for the Cholesky decompositions? Default is TRUE
#' @param ... Currently not used.
#' @return A cov_rSPDE object
#' @method simulate cov_rSPDEobj
#' @export 

simulate.cov_rSPDEobj <- function(object, nsim = 1,
                                user_nu = NULL,
                                user_kappa = NULL,
                                user_tau = NULL,
                                user_m = NULL,
                                pivot = TRUE,
                                ...)
{
  
  object <- update.cov_rSPDEobj(object,
                                user_nu,
                                user_kappa,
                                user_tau,
                                user_m)

  d <- object$d
  kappa <- object$kappa
  nu <- object$nu
  tau <- object$tau
  m <- object$m
  
  alpha <- nu + d/2
  m_alpha <- max(1,floor(alpha))
  
  fem_mesh_matrices <- object$fem_mesh_matrices

  
  ## simulation 
  if (alpha%%1==0){# simulation in integer case
    Q <- rspde_Prec_int(kappa, nu, tau, d, fem_mesh_matrices)
    Z = rnorm(sizeL * nsim)
    dim(Z) <- c(sizeL, nsim)
    if(pivot){
      LQ = chol(Q, pivot = TRUE)
      reorder = attr(LQ,"pivot")
      X = solve(LQ,Z)
      # order back
      orderback = numeric(length(reorder))
      orderback[reorder] = 1:length(reorder)
      X = X[orderback]
    } else{
      LQ = chol(Q)
      X = solve(LQ,Z)
    }
  }else{# simulation in non-integer case
    ## get rational approximation coefficients: r,p,k
    # name of the table we need to reach
    
    L = fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/kappa^2
    sizeL = dim(L)[1]
    
    mt <- get(paste0("m", m, "t"))
    # get r
    r = sapply(1:(m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    # get p
    p = sapply(1:(m), function(i) {
      approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
    })
    # get k
    k = approx(mt$nu, mt$k, cut_decimals(nu))$y
    # one part to construct Q_i, i<=m
    if (m_alpha == 1){
      Malpha =  (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
    } else if (m_alpha > 1){
      Malpha =  fem_mesh_matrices[["c0"]] + m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
      for(j in 2:m_alpha){
        Malpha = Malpha +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
    
    # another part to construct Q_i, i<=m
    if (m_alpha == 1){
      Malpha2 =  (fem_mesh_matrices[["g1"]] + fem_mesh_matrices[["g2"]]/(kappa^2))
    } else if (m_alpha > 1){
      Malpha2 =  fem_mesh_matrices[["g1"]] + m_alpha * fem_mesh_matrices[["g2"]]/(kappa^2)
      for(j in 2:m_alpha){
        Malpha2 = Malpha2 +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j+1)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
    
    # initialize the simulated uhat as X = 0
    X = 0
    # compute sum of x_i, i from 1 to m
    for (i in 1:m){
      # one sample from standard  normal multivariate gaussian
      Z = rnorm(sizeL * nsim)
      dim(Z) <- c(sizeL, nsim)
      # cholesky decomposition of Q_i
      # R is an upper triangle matrix
      if(pivot){
        R = chol(1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha), pivot = TRUE)
        reorder = attr(R,"pivot")
        # order back
        orderback = numeric(length(reorder))
        orderback[reorder] = 1:length(reorder)
        X = X + solve(R,Z)[orderback]
      } else{
        R = chol(1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
        # use backsolve for upper triangle matrix
        X = X + solve(R,Z)
      }
    }
    # now let us add the k_part which is Q_m+1
    Z = rnorm(sizeL * nsim)
    dim(Z) <- c(sizeL, nsim)
    count = m_alpha
    diag_C0 <- diag(fem_mesh_matrices[["c0"]])
    
    if (m_alpha%%2 == 0){# when m_alpha is even
      CR = sqrt(k*diag_C0)
      k_part = CR * Z
      while (count > 0){
        k_part = solve(L, diag_C0 * k_part)
        count = count-2
      }
    }else{
      if(pivot){
        LR <- chol((1/k)*L, pivot=TRUE)
        reorder <- attr(LR,"pivot")
        orderback = numeric(length(reorder))
        orderback[reorder] = 1:length(reorder)
        k_part = solve(LR,Z)[orderback]
      } else{
        LR = chol((1/k)*L)
        k_part = solve(LR,Z)
      }

      while (count > 1){
        k_part = solve(L,diag_C0 * k_part)
        count = count-2
      }
    }
    X = X + k_part
    # compute u, the simulated field
    X = kappa^(-alpha)*X
    X = X/tau
  }
  return(X)
}


#' Prediction of a fractional SPDE using a rational SPDE approximation
#'
#' The function is used for computing kriging predictions based on data \eqn{Y_i = u(s_i) + \epsilon_i},
#' where \eqn{\epsilon}{\epsilon} is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by
#' a fractional SPDE \eqn{L^\beta u(s) = W}{L^\beta u(s) = W}, where \eqn{W}{W} is Gaussian white noise.
#'
#' @param object The rational SPDE approximation, computed using \code{\link{fractional.operators}},
#' \code{\link{matern.operators}}, or \code{\link{spde.matern.operators}}.
#' @param A A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
#' @param Aprd A matrix linking the prediction locations to the basis of the FEM approximation of the latent model.
#' @param Y A vector with the observed data, can also be a matrix where the columns are observations
#' of independent replicates of \eqn{u}.
#' @param sigma.e The standard deviation of the Gaussian measurement noise. Put to zero if the model
#' does not have measurement noise.
#' @param compute.variances Set to also TRUE to compute the kriging variances.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict rSPDEobj
#'
#' @examples
#' #Sample a Gaussian Matern process on R using a rational approximation
#' kappa <- 10
#' sigma <- 1
#' nu <- 0.8
#' sigma.e <- 0.3
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x <- seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation
#' op <- matern.operators(kappa = kappa, sigma = sigma,
#'                        nu = nu, G=fem$G, C = fem$C, d = 1)
#'
#' #Sample the model
#' u <- simulate(op)
#'
#' #Create some data
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y <- as.vector(A%*%u + sigma.e*rnorm(10))
#'
#' #compute kriging predictions at the FEM grid
#' A.krig <- rSPDE.A1d(x, x)
#' u.krig <- predict(op, A = A, Aprd = A.krig, Y = Y, sigma.e = sigma.e,
#'                   compute.variances= TRUE)
#'
#' plot(obs.loc, Y, ylab = "u(x)", xlab = "x", main = "Data and prediction",
#'      ylim = c(min(u.krig$mean - 2*sqrt(u.krig$variance)),
#'               max(u.krig$mean + 2*sqrt(u.krig$variance))))
#' lines(x, u.krig$mean)
#' lines(x, u.krig$mean + 2*sqrt(u.krig$variance), col = 2)
#' lines(x, u.krig$mean - 2*sqrt(u.krig$variance), col = 2)

predict.rSPDEobj <- function(object, A, Aprd, Y, sigma.e, compute.variances = FALSE,...)
{
  Y <- as.matrix(Y)
  if (dim(Y)[1] != dim(A)[1])
    stop("the dimensions of A does not match the number of observations")
  
  out <- list()
  if(length(sigma.e) == 1){
    if (sigma.e < 0) {
      stop("sigma.e must be non-negative")
    } else if (sigma.e > 0) {
      A <- A %*% object$Pr
      AA <- Aprd %*% object$Pr
      Qhat <- object$Q + (t(A) %*% A) / sigma.e^2
      out$mean <- as.matrix(AA %*% solve(Qhat, t(A) %*% Y / sigma.e^2))
      if (compute.variances) {
        out$variance <- diag(AA %*% solve(Qhat,t(AA)))
      }
    } else { #no nugget
      Ahat <- A %*% object$Pr
      QiAt <- solve(object$Q, t(Ahat))
      AQiA <- Ahat %*% QiAt
      xhat <- solve(object$Q, t(Ahat)%*%solve(AQiA, Y))
      out$mean <- as.vector(Aprd %*% xhat)
      if (compute.variances) {
        AA <- Aprd %*% object$Pr
        M <- object$Q - QiAt %*% solve(AQiA, t(QiAt))
        out$variance <- diag(AA%*%M%*%t(AA))
      }
    }  
  } else if(dim(Y)[1] == length(sigma.e)){
    Q.e <- Diagonal(length(sigma.e),1/sigma.e^2)
    A <- A %*% object$Pr
    AA <- Aprd %*% object$Pr
    Qhat <- object$Q + t(A) %*% Q.e%*% A 
    out$mean <- as.matrix(AA %*% solve(Qhat, t(A) %*%Q.e%*% Y))
    if (compute.variances) {
      out$variance <- diag(AA %*% solve(Qhat,t(AA)))
    }
  }
  return(out)
}

#' Log-likelihood function for latent Gaussian fractional SPDE model
#'
#' This function evaluates the log-likelihood function for a fractional SPDE model
#' \eqn{L^\beta u(s) = W}{L^\beta u(s) = W} that is observed under Gaussian measurement
#' noise: \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = x(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i}
#' are iid mean-zero Gaussian variables and \eqn{x(s) = \mu(s) + u(s)}{x(s) = \mu(s) + u(s)}, where 
#' \eqn{\mu(s)}{\mu(s)} is the expectation vector of the latent field.
#'
#' @param obj The rational SPDE approximation, computed using \code{\link{fractional.operators}},
#' \code{\link{matern.operators}}, or \code{\link{spde.matern.operators}}.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0). 
#' @return The log-likelihood value.
#' @export
#' @note This example below shows how the function can be used to evaluate the likelihood of a latent
#' Matern model. Se \code{\link{matern.loglike}} for an example of how this can be used for maximum
#' likelihood estimation.
#' @seealso \code{\link{matern.loglike}}, \code{\link{spde.matern.loglike}}
#'
#' @examples
#' #Sample a Gaussian Matern process on R using a rational approximation
#' kappa = 10
#' sigma = 1
#' nu = 0.8
#' sigma.e = 0.3
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation
#' op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
#'                        G = fem$G, C = fem$C, d = 1)
#'
#' #Sample the model
#' u <- simulate(op)
#'
#' #Create some data
#' obs.loc <- runif(n = 10, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' Y = as.vector(A%*%u + sigma.e*rnorm(10))
#'
#' #compute log-likelihood of the data
#' lik1 <- rSPDE.loglike(op, Y, A, sigma.e)
#' cat(lik1)

rSPDE.loglike <- function(obj, Y, A, sigma.e, mu=0)
{
  Y = as.matrix(Y)
  if (length(dim(Y)) == 2) {
    n.rep = dim(Y)[2]
    n = dim(Y)[1]
  } else {
    n.rep = 1
    if (length(dim(Y)) == 1){
      n = dim(Y)[1]
    } else {
      n = length(Y)
    }
  }
  if(length(sigma.e)==1){
    Q.e <- Diagonal(n)/sigma.e^2
    nugget = rep(sigma.e^2,n)
  } else {
    if(length(sigma.e) != n){
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e),1/sigma.e^2)
    nugget = sigma.e^2
  }
  R = Matrix::Cholesky(obj$Pl)
  prior.ld = 4 * c(determinant(R, logarithm = TRUE)$modulus) - sum(log(diag(obj$C)))
  
  A = A %*% obj$Pr
  Q.post = obj$Q + t(A) %*% Q.e %*% A 
  R.post = Matrix::Cholesky(Q.post)
  posterior.ld = 2 * c(determinant(R.post, logarithm = TRUE)$modulus)
  
  AtY = t(A) %*% Q.e %*% Y 
  mu.post <- mu + solve(R.post, AtY, system = "A")
  
  lik = n.rep * (prior.ld - posterior.ld - dim(A)[1] * log(2*pi) - sum(log(nugget))) / 2
  
  if (n.rep > 1) {
    lik = lik - 0.5 * sum(colSums((mu.post -mu)* (obj$Q %*% (mu.post - mu))))
    v = Q.e%*%(Y - A %*% mu.post)
    lik = lik - 0.5 * sum(colSums((Y - A %*% mu.post)*v)) 
  } else {
    lik = lik - 0.5 * (t(mu.post - mu) %*% obj$Q %*% (mu.post - mu) + t(Y - A %*% mu.post) %*% Q.e %*% (Y - A %*% mu.post))
  }
  return(as.double(lik))
}

#' @name cov_rSPDE_matern.loglike
#' @title Log-likelihood function for latent Gaussian fractional SPDE model
#' @description This function evaluates the log-likelihood function for covariance-based
#' rational approximation of a fractional SPDE model
#' @param object A cov_rSPDE object
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location to the finite element basis.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0). 
#' @param user_kappa If non-null, update the range parameter of the covariance function.
#' @param user_tau If non-null, update the standard deviation of the covariance function.
#' @param user_nu If non-null, update the shape parameter of the covariance function.
#' @param user_m If non-null, update the order of the rational approximation, which needs to be a positive integer.
#' @param pivot Should pivoting be used for the Cholesky decompositions? Default is TRUE
#' @return A cov_rSPDE object
#' @export

cov_rSPDE_matern.loglike <- function(object, Y, A, sigma.e, mu=0,
                            user_nu = NULL,
                            user_kappa = NULL,
                            user_tau = NULL,
                            user_m = NULL,
                            pivot=TRUE)
{
 
  Y = as.matrix(Y)
  if (length(dim(Y)) == 2) {
    n.rep = dim(Y)[2]
    n = dim(Y)[1]
  } else {
    n.rep = 1
    if (length(dim(Y)) == 1){
      n = dim(Y)[1]
    } else {
      n = length(Y)
    }
  }
  
  ## get relevant parameters
  object <- update.cov_rSPDEobj(object,
                                user_nu,
                                user_kappa,
                                user_tau,
                                user_m)
  
  d <- object$d
  kappa <- object$kappa
  nu <- object$nu
  tau <- object$tau
  m <- object$m
  
  fem_mesh_matrices <- object$fem_mesh_matrices
  
  if(length(sigma.e)==1){
    Q.e <- Diagonal(n)/sigma.e^2
    nugget = rep(sigma.e^2,n)
  } else {
    if(length(sigma.e) != n){
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e),1/sigma.e^2)
    nugget = sigma.e^2
  }
  
  alpha <- nu + d/2
  m_alpha = max(1,floor(alpha))
  ## get L matrix
  L = fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/kappa^2
  sizeL = dim(L)[1]
  
  if (alpha%%1==0){# loglikelihood in integer case
    ## construct Q
    Q <- rspde_Prec_int(kappa, nu, tau, d, fem_mesh_matrices)
    
    R = chol(Q,pivot = pivot)
    logQ = 2*sum(log(diag(R)))
    ## compute Q_x|y
    Q_xgiveny = (t(A)%*% Q.e %*% A) + Q
    ## construct mu_x|y 
    mu_xgiveny = t(A) %*% Q.e %*% Y
    # upper triangle with reordering
    R = chol(Q_xgiveny, pivot = pivot)
    if(pivot){
      reorder = attr(R,"pivot")
      # make it lower triangle
      R = t(R)
      v = solve(R,mu_xgiveny[reorder])
      mu_xgiveny = solve(t(R),v)
      # order back
      orderback = numeric(length(reorder))
      orderback[reorder] = 1:length(reorder)
      mu_xgiveny = mu_xgiveny[orderback]
    } else{
      R = t(R)
      v = solve(R,mu_xgiveny)
      mu_xgiveny = solve(t(R),v)
    }
    
    mu_xgiveny <- mu + mu_xgiveny

    ## compute log|Q_xgiveny|
    log_Q_xgiveny = 2*sum(log(diag(R)))
    ## compute mu_x|y*Q*mu_x|y
    if(n.rep>1){
      mu_part = sum(colSums((mu_xgiveny -mu) * (Q %*% (mu_xgiveny-mu))))
    } else{
      mu_part = t(mu_xgiveny -mu) %*% Q %*% (mu_xgiveny-mu)      
    }

    ## compute central part
    if(n.rep>1){
      central_part = sum(colSums((Y-A %*% mu_xgiveny) *(Q.e %*% (Y-A %*% mu_xgiveny))))
    } else{
      central_part = t(Y-A %*% mu_xgiveny) %*% Q.e %*% (Y-A %*% mu_xgiveny) 
    }
    
    ## compute log|Q_epsilon|
    log_Q_epsilon = -sum(log(nugget))


    ## wrap up
    log_likelihood = n.rep * (logQ + log_Q_epsilon - log_Q_xgiveny) - mu_part - central_part
    if(n.rep>1){
      log_likelihood = log_likelihood - dim(A)[1] * n.rep*log(2*pi)
    } else{
      log_likelihood = log_likelihood - length(Y)*log(2*pi)
    }

    log_likelihood = log_likelihood/2
    
  }else{# loglikelihood in non-integer case
    ## get rational approximation coefficients: r,p,k
    # name of the table we need to reach
    mt <- get(paste0("m", m, "t"))
    # get r
    r = sapply(1:(m), function(i) {
      approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
    })
    # get p
    p = sapply(1:(m), function(i) {
      approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
    })
    # get k
    k = approx(mt$nu, mt$k, cut_decimals(nu))$y
    # one part to construct Q_i, i<=m
    if (m_alpha == 1){
      Malpha =  (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
    } else if (m_alpha > 1){
      Malpha =  fem_mesh_matrices[["c0"]] + m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
      for(j in 2:m_alpha){
        Malpha = Malpha +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
    
    # another part to construct Q_i, i<=m
    if (m_alpha == 1){
      Malpha2 =  (fem_mesh_matrices[["g1"]] + fem_mesh_matrices[["g2"]]/(kappa^2))
    } else if (m_alpha > 1){
      Malpha2 =  fem_mesh_matrices[["g1"]] + m_alpha * fem_mesh_matrices[["g2"]]/(kappa^2)
      for(j in 2:m_alpha){
        Malpha2 = Malpha2 +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j+1)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
    
    ## compute log|Q|
    logQ = 0
    # compute sum of x_i, i from 1 to m
    for (i in 1:m){
      
      # cholesky decomposition of Q_i
      # R is an upper triangle matrix
      R = chol(1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
      logQ = logQ + sum(log(diag(R)))
    }
    # now let us add the logk_part which is log|Q_m+1|
    logkpart = 0
    LR = chol(L)
    loglpart = sum(log(diag(LR)))
    logcinvpart = -sum(log(diag(fem_mesh_matrices[["c0"]])))
    logkpart = m_alpha*loglpart + (m_alpha-1)*logcinvpart -log(k)
    # get logQ
    logQ = 2*logQ + logkpart
    logQ = logQ + (m+1)*sizeL*log(kappa ^ (2 * alpha) + tau ^ 2)
    
    ## construct Q
    Q = 1/r[1] * (Malpha + Malpha2/kappa^2 - p[1] * Malpha)
    
    if(length(r)>1){
      for (i in 2:length(r)) {
        Q = bdiag(Q, 1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
      }
    }
    
    # add k_part into Q
    
    
    if (m_alpha == 1){
      Kpart = 1/k * (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
    } else if (m_alpha > 1){
      Kpart = 1/k * fem_mesh_matrices[["c0"]] + 1/k * m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
      for(j in 2:m_alpha){
        Kpart = Kpart + 1/k * choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
      }
    } else {
      stop("Something is wrong with the value of nu!")
    }
    
    Q = bdiag(Q, Kpart)
    
    Q = Q * kappa ^ (4 * beta)
    
    Q = tau ^ 2 * Q
    
    ## compute Q_x|y
    Q_xgiveny = kronecker(matrix(1,length(r)+1,length(r)+1),t(A)%*% Q.e %*% A) + Q
    ## construct mu_x|y 
    Abar = kronecker(matrix(1,1,length(r)+1),A)
    mu_xgiveny = t(Abar) %*% Q.e %*% Y
    # upper triangle with reordering
    R = chol(Q_xgiveny, pivot = pivot)
    if(pivot){
      reorder = attr(R,"pivot")
      # make it lower triangle
      R = t(R)
      v = solve(R,mu_xgiveny[reorder])
      mu_xgiveny = solve(t(R),v)
      
      orderback = numeric(length(reorder))
      orderback[reorder] = 1:length(reorder)
      mu_xgiveny = mu_xgiveny[orderback]
    } else{
      R = t(R)
      v = solve(R,mu_xgiveny)
      mu_xgiveny = solve(t(R),v)
    }
    mu_xgiveny <- mu + mu_xgiveny

    ## compute log|Q_xgiveny|
    log_Q_xgiveny = 2*sum(log(diag(R)))
    ## compute mu_x|y*Q*mu_x|y
    if(n.rep>1){
      mu_part = sum(colSums((mu_xgiveny -mu) * (Q %*% (mu_xgiveny-mu))))
    } else{
      mu_part = t(mu_xgiveny -mu) %*% Q %*% (mu_xgiveny-mu)      
    }
    ## compute central part
    if(n.rep>1){
      central_part = sum(colSums((Y-Abar %*% mu_xgiveny) *(Q.e %*% (Y-Abar %*% mu_xgiveny))))
    } else{
      central_part = t(Y-A %*% mu_xgiveny) %*% Q.e %*% (Y-A %*% mu_xgiveny) 
    }
    ## compute log|Q_epsilon|
    log_Q_epsilon = -sum(log(nugget))
    ## wrap up
    log_likelihood = n.rep * (logQ + log_Q_epsilon - log_Q_xgiveny) - mu_part - central_part
    if(n.rep>1){
      log_likelihood = log_likelihood - dim(A)[1] * n.rep*log(2*pi)
    } else{
      log_likelihood = log_likelihood - length(Y)*log(2*pi)
    }
    log_likelihood = log_likelihood/2
  }
  
  return(log_likelihood)
  
}





#' Log-likelihood for a latent Gaussian Matern model using a rational SPDE approximation
#'
#' This function evaluates the log-likelihood function for a Gaussian process with a Matern covariance
#' function, that is observed under Gaussian measurement noise:
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i} are
#' iid mean-zero Gaussian variables. The latent model is approximated using a rational approximation
#' of the fractional SPDE model corresponding to the Gaussian process.
#'
#' @param kappa Range parameter of the latent process.
#' @param sigma Standard deviation of the latent process.
#' @param nu Shape parameter of the latent process.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization of the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a positive integer.
#' The default value is 1.
#'
#' @return The log-likelihood value.
#' @export
#' @seealso \code{\link{spde.matern.loglike}}, \code{\link{rSPDE.loglike}}, \code{\link{matern.operators}}.
#'
#' @examples
#' #this example illustrates how the function can be used for maximum likelihood estimation
#' set.seed(123)
#' #Sample a Gaussian Matern process on R using a rational approximation
#' sigma = 1
#' nu = 0.8
#' kappa = 1
#' sigma.e = 0.3
#' n.rep = 10
#' n.obs = 100
#' n.x = 51
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation
#' op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
#'                        G = fem$G, C = fem$C, d = 1)
#'
#' #Sample the model
#' u <- simulate(op, n.rep)
#'
#' #Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs*n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y = as.matrix(A%*%u + sigma.e*noise)
#'
#' #define negative likelihood function for optimization using matern.loglike
#' mlik <- function(theta, Y, G, C, A){
#' return(-matern.loglike(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]),
#'                        Y = Y, G = G, C = C, A = A, d = 1))
#' }
#'
#' #The parameters can now be estimated by maximizing mlik with optim
#' \donttest{
#' #Choose some reasonable starting values depending on the size of the domain
#' theta0 = log(c(sqrt(8), sqrt(var(c(Y))), 0.9, 0.01))
#'
#' #run estimation and display the results
#' theta <- optim(theta0, mlik, Y = Y, G = fem$G, C = fem$C, A = A)
#'
#' print(data.frame(kappa = c(kappa,exp(theta$par[1])), sigma = c(sigma,exp(theta$par[2])),
#'                  nu = c(nu,exp(theta$par[3])), sigma.e = c(sigma.e,exp(theta$par[4])),
#'                  row.names = c("Truth","Estimates")))
#' }

matern.loglike <- function(kappa,
                           sigma,
                           nu,
                           sigma.e,
                           Y,
                           G,
                           C,
                           A,
                           d = 2,
                           m = 1)
{
  op <- matern.operators(kappa, sigma, nu, G, C, d = d, m = m)
  return(rSPDE.loglike(op, Y, A, sigma.e))
}


#' @name cov_rSPDE_matern.loglike
#' @title Log-likelihood function for latent Gaussian fractional SPDE model
#' @description This function evaluates the log-likelihood function for covariance-based
#' rational approximation of a fractional SPDE model
#' @param kappa Range parameter of the latent process.
#' @param tau Standard deviation of the latent process.
#' @param nu Shape parameter of the latent process.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param mu Expectation vector of the latent field (default = 0). 
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization of the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a positive integer.
#' The default value is 2.
#' @param pivot Should pivoting be used for the Cholesky decompositions? Default is TRUE
#' @return A cov_rSPDE object
#' @export

matern.cov_rSPDE.loglike <- function(kappa,
                                     tau,
                                     nu,
                                     sigma.e,
                                     mu=0,
                                     Y,
                                     G,
                                     C,
                                     A,
                                     d = 2,
                                     m = 2,
                                     pivot=TRUE)
{
  obj_cov_rSPDE <- matern.cov_rSPDE.operators(C=C,
                                                          G=G,
                                                          nu=nu,
                                                          kappa=kappa,
                                                          tau=tau,
                                                          m=m,
                                                          d=d)
  
  return(cov_rSPDE_matern.loglike(object=obj_cov_rSPDE, Y=Y, A=A, sigma.e=sigma.e, mu=mu,
                       user_nu = NULL,
                       user_kappa = NULL,
                       user_tau = NULL,
                       user_m = NULL,
                       pivot=pivot))
}






#' Log-likelihood for a latent Gaussian Matern SPDE model using a rational SPDE approximation
#'
#' This function evaluates the log-likelihood function for observations of a Gaussian process defined as
#' the solution to the SPDE \deqn{(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W},
#'
#' The observations are assumed to be generated as
#' \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i} are
#' iid mean-zero Gaussian variables. The latent model is approximated using a rational approximation
#' of the fractional SPDE model.
#'
#' @param kappa Vector with the, possibly spatially varying, range parameter evaluated at the locations
#' of the mesh used for the finite element discretization of the SPDE.
#' @param tau Vector with the, possibly spatially varying, precision parameter evaluated at the locations
#' of the mesh used for the finite element discretization of the SPDE.
#' @param nu Shape parameter of the covariance function, related to \eqn{\beta} through the equation
#' \eqn{\beta = (\nu + d/2)/2}.
#' @param sigma.e The standard deviation of the measurement noise.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param G The stiffness matrix of a finite element discretization of the domain.
#' @param C The mass matrix of a finite element discretization of the domain.
#' @param A A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a positive integer.
#' The default value is 1.
#'
#' @return The log-likelihood value.
#' @export
#' @seealso \code{\link{matern.loglike}}, \code{\link{rSPDE.loglike}}.
#'
#' @examples
#' #this example illustrates how the function can be used for maximum likelihood estimation
#' set.seed(1)
#' #Sample a Gaussian Matern process on R using a rational approximation
#' sigma.e = 0.1
#' n.rep = 10
#' n.obs = 100
#' n.x = 51
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = n.x)
#' fem <- rSPDE.fem1d(x)
#'
#' tau = rep(0.5,n.x)
#' nu = 0.8
#' kappa = rep(1,n.x)
#'
#' #compute rational approximation
#' op <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu,
#'                             G = fem$G, C = fem$C, d = 1)
#'
#' #Sample the model
#' u <- simulate(op, n.rep)
#'
#' #Create some data
#' obs.loc <- runif(n = n.obs, min = 0, max = 1)
#' A <- rSPDE.A1d(x, obs.loc)
#' noise <- rnorm(n.obs*n.rep)
#' dim(noise) <- c(n.obs, n.rep)
#' Y = as.matrix(A%*%u + sigma.e*noise)
#'
#' #define negative likelihood function for optimization using matern.loglike
#' mlik <- function(theta, Y, G, C, A){
#' return(-spde.matern.loglike(rep(exp(theta[1]),n.x), rep(exp(theta[2]),n.x),
#'                             exp(theta[3]), exp(theta[4]),
#'                             Y = Y, G = G, C = C, A = A, d = 1))
#' }
#'
#'#' #The parameters can now be estimated by maximizing mlik with optim
#' \donttest{
#' #Choose some reasonable starting values depending on the size of the domain
#' theta0 = log(c(sqrt(8), 1/sqrt(var(c(Y))), 0.9, 0.01))
#'
#' #run estimation and display the results
#' theta <- optim(theta0, mlik, Y = Y, G = fem$G, C = fem$C, A = A)
#'
#' print(data.frame(kappa = c(kappa[1],exp(theta$par[1])), tau = c(tau[1],exp(theta$par[2])),
#'                  nu = c(nu,exp(theta$par[3])), sigma.e = c(sigma.e,exp(theta$par[4])),
#'                  row.names = c("Truth","Estimates")))
#'}

spde.matern.loglike <- function(kappa,
                                tau,
                                nu,
                                sigma.e,
                                Y,
                                G,
                                C,
                                A,
                                d = 2,
                                m = 1)
{
  op <- spde.matern.operators(kappa, tau, nu, G, C, d = d, m = m)
  return(rSPDE.loglike(op, Y, A, sigma.e))
}



#' @name predict.cov_rSPDEobj
#' @title Prediction of a fractional SPDE using the covariance-based rational SPDE approximation
#' @description The function is used for computing kriging predictions based on data \eqn{Y_i = u(s_i) + \epsilon_i},
#' where \eqn{\epsilon}{\epsilon} is mean-zero Gaussian measurement noise and \eqn{u(s)}{u(s)} is defined by
#' a fractional SPDE \eqn{L^\beta u(s) = W}{L^\beta u(s) = W}, where \eqn{W}{W} is Gaussian white noise.
#' @param object A cov_rSPDE object
#' @param A A matrix linking the measurement locations to the basis of the FEM approximation of the latent model.
#' @param Aprd A matrix linking the prediction locations to the basis of the FEM approximation of the latent model.
#' @param Y A vector with the observed data, can also be a matrix where the columns are observations
#' of independent replicates of \eqn{u}.
#' @param sigma.e The standard deviation of the Gaussian measurement noise. Put to zero if the model
#' does not have measurement noise.
#' @param mu Expectation vector of the latent field (default = 0). 
#' @param compute.variances Set to also TRUE to compute the kriging variances.
#' @param pivot Should pivoting be used on the Cholesky decompositions?
#' @param ... further arguments passed to or from other methods.
#' @return A list with elements
#' \item{mean }{The kriging predictor (the posterior mean of u|Y).}
#' \item{variance }{The posterior variances (if computed).}
#' @export
#' @method predict cov_rSPDEobj


predict.cov_rSPDEobj <- function(object, A, Aprd, Y, sigma.e, mu=0, 
                                 compute.variances = FALSE, pivot=TRUE,
                                 ...)
{
  Y <- as.matrix(Y)
  if (dim(Y)[1] != dim(A)[1])
    stop("the dimensions of A does not match the number of observations")
  
  out <- list()
  
  d <- object$d
  kappa <- object$kappa
  nu <- object$nu
  tau <- object$tau
  m <- object$m
  
  fem_mesh_matrices <- object$fem_mesh_matrices
  
  no_nugget = FALSE
  
  if(length(sigma.e)==1){
    if(sigma.e==0){
      no_nugget = TRUE
    } else{
      Q.e <- Diagonal(n)/sigma.e^2
      nugget = rep(sigma.e^2,n)
    }
  } else {
    if(length(sigma.e) != n){
      stop("the length of sigma.e does not match the number of observations")
    }
    Q.e <- Diagonal(length(sigma.e),1/sigma.e^2)
    nugget = sigma.e^2
  }
  
  alpha <- nu + d/2
  m_alpha = max(1,floor(alpha))
  ## get L matrix
  L = fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/kappa^2
  sizeL = dim(L)[1]
  
  if(!no_nugget){
    if (alpha%%1==0){# loglikelihood in integer case
      ## construct Q
      Q <- rspde_Prec_int(kappa, nu, tau, d, fem_mesh_matrices)
      
      R = chol(Q,pivot = pivot)
      logQ = 2*sum(log(diag(R)))
      ## compute Q_x|y
      Q_xgiveny = (t(A)%*% Q.e %*% A) + Q
      ## construct mu_x|y 
      mu_xgiveny = t(A) %*% Q.e %*% Y
      # upper triangle with reordering
      R = chol(Q_xgiveny, pivot = pivot)
      if(pivot){
        reorder = attr(R,"pivot")
        # make it lower triangle
        R = t(R)
        v = solve(R,mu_xgiveny[reorder])
        mu_xgiveny = solve(t(R),v)
        # order back
        orderback = numeric(length(reorder))
        orderback[reorder] = 1:length(reorder)
        mu_xgiveny = mu_xgiveny[orderback]
      } else{
        R = t(R)
        v = solve(R,mu_xgiveny)
        mu_xgiveny = solve(t(R),v)
      }
      
      mu_xgiveny <- mu + mu_xgiveny
      out$mean <- Aprd %*% mu_xgiveny

      if (compute.variances) {
        out$variance <- diag(Aprd %*% solve(Q_xgiveny,t(Aprd)))
      }
      
    }else{# loglikelihood in non-integer case
      ## get rational approximation coefficients: r,p,k
      # name of the table we need to reach
      mt <- get(paste0("m", m, "t"))
      # get r
      r = sapply(1:(m), function(i) {
        approx(mt$nu, mt[[paste0("r", i)]], cut_decimals(nu))$y
      })
      # get p
      p = sapply(1:(m), function(i) {
        approx(mt$nu, mt[[paste0("p", i)]], cut_decimals(nu))$y
      })
      # get k
      k = approx(mt$nu, mt$k, cut_decimals(nu))$y
      # one part to construct Q_i, i<=m
      if (m_alpha == 1){
        Malpha =  (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
      } else if (m_alpha > 1){
        Malpha =  fem_mesh_matrices[["c0"]] + m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
        for(j in 2:m_alpha){
          Malpha = Malpha +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
      
      # another part to construct Q_i, i<=m
      if (m_alpha == 1){
        Malpha2 =  (fem_mesh_matrices[["g1"]] + fem_mesh_matrices[["g2"]]/(kappa^2))
      } else if (m_alpha > 1){
        Malpha2 =  fem_mesh_matrices[["g1"]] + m_alpha * fem_mesh_matrices[["g2"]]/(kappa^2)
        for(j in 2:m_alpha){
          Malpha2 = Malpha2 +  choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j+1)]]/(kappa^(2*j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
      
      ## compute log|Q|
      logQ = 0
      # compute sum of x_i, i from 1 to m
      for (i in 1:m){
        
        # cholesky decomposition of Q_i
        # R is an upper triangle matrix
        R = chol(1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
        logQ = logQ + sum(log(diag(R)))
      }
      # now let us add the logk_part which is log|Q_m+1|
      logkpart = 0
      LR = chol(L)
      loglpart = sum(log(diag(LR)))
      logcinvpart = -sum(log(diag(fem_mesh_matrices[["c0"]])))
      logkpart = m_alpha*loglpart + (m_alpha-1)*logcinvpart -log(k)
      # get logQ
      logQ = 2*logQ + logkpart
      logQ = logQ + (m+1)*sizeL*log(kappa ^ (2 * alpha) + tau ^ 2)
      
      ## construct Q
      Q = 1/r[1] * (Malpha + Malpha2/kappa^2 - p[1] * Malpha)
      
      if(length(r)>1){
        for (i in 2:length(r)) {
          Q = bdiag(Q, 1/r[i] * (Malpha + Malpha2/kappa^2 - p[i] * Malpha))
        }
      }
      
      # add k_part into Q
      
      
      if (m_alpha == 1){
        Kpart = 1/k * (fem_mesh_matrices[["c0"]] + fem_mesh_matrices[["g1"]]/(kappa^2))
      } else if (m_alpha > 1){
        Kpart = 1/k * fem_mesh_matrices[["c0"]] + 1/k * m_alpha * fem_mesh_matrices[["g1"]]/(kappa^2)
        for(j in 2:m_alpha){
          Kpart = Kpart + 1/k * choose(m_alpha, j) * fem_mesh_matrices[[paste0("g",j)]]/(kappa^(2*j))
        }
      } else {
        stop("Something is wrong with the value of nu!")
      }
      
      Q = bdiag(Q, Kpart)
      
      Q = Q * kappa ^ (4 * beta)
      
      Q = tau ^ 2 * Q
      
      ## compute Q_x|y
      Q_xgiveny = kronecker(matrix(1,length(r)+1,length(r)+1),t(A)%*% Q.e %*% A) + Q
      ## construct mu_x|y 
      Abar = kronecker(matrix(1,1,length(r)+1),A)
      mu_xgiveny = t(Abar) %*% Q.e %*% Y
      # upper triangle with reordering
      R = chol(Q_xgiveny, pivot = pivot)
      if(pivot){
        reorder = attr(R,"pivot")
        # make it lower triangle
        R = t(R)
        v = solve(R,mu_xgiveny[reorder])
        mu_xgiveny = solve(t(R),v)
        
        orderback = numeric(length(reorder))
        orderback[reorder] = 1:length(reorder)
        mu_xgiveny = mu_xgiveny[orderback]
      } else{
        R = t(R)
        v = solve(R,mu_xgiveny)
        mu_xgiveny = solve(t(R),v)
      }
      mu_xgiveny <- mu + mu_xgiveny
      
      Aprd_bar = kronecker(matrix(1,1,length(r)+1),Aprd)
      
      out$mean <- Aprd_bar %*% mu_xgiveny
      
      if (compute.variances) {
        out$variance <- diag(Aprd_bar %*% solve(Q_xgiveny,t(Aprd_bar)))
      }
    }
  } else{
    integer_alpha <- (alpha%%1 == 0)
    if(integer_alpha){
      Abar = A
      Aprd_bar = Aprd
      Q <- rspde_Prec_int(kappa=kappa, nu=nu, tau=tau, d=d, 
                          fem_mesh_matrices=fem_mesh_matrices)
    } else{
      Abar = kronecker(matrix(1,1,rspde_order+1),A)
      Aprd_bar = kronecker(matrix(1,1,rspde_order+1),Aprd)
      Q <- rspde_Prec(kappa=kappa, nu=nu, tau=tau, rspde_order=m, d=d, 
                      fem_mesh_matrices=fem_mesh_matrices)
    }
  
    QiAt <- solve(Q, t(Abar))
    AQiA <- Abar %*% QiAt
    xhat <- solve(Q, t(Abar)%*%solve(AQiA, Y))
    
    out$mean <- as.vector(Aprd_bar %*% xhat)
    if (compute.variances) {
      M <- Q - QiAt %*% solve(AQiA, t(QiAt))
      out$variance <- diag(Aprd_bar%*%M%*%t(Aprd_bar))
    }
  }
  
 
  return(out)
}
