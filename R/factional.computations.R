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
simulate <- function(object, nsim) {
  UseMethod("simulate", object)
}

#' Simulation of a fractional SPDE using a rational SPDE approximation
#'
#' The function samples a Gaussian random field based on a pre-computed rational SPDE approximation.
#'
#' @param object The rational SPDE approximation, computed using \code{\link{fractional.operators}},
#' \code{\link{matern.operators}}, or \code{\link{spde.matern.operators}}.
#' @param nsim The number of simulations.
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

simulate.rSPDEobj <- function(object, nsim = 1)
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
#' Y <- as.vector(A\%*\%u + sigma.e*rnorm(10))
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
#' noise: \eqn{Y_i = u(s_i) + \epsilon_i}{Y_i = u(s_i) + \epsilon_i}, where \eqn{\epsilon_i}{\epsilon_i}
#' are iid mean-zero Gaussian variables.
#'
#' @param obj The rational SPDE approximation, computed using \code{\link{fractional.operators}},
#' \code{\link{matern.operators}}, or \code{\link{spde.matern.operators}}.
#' @param Y The observations, either a vector or a matrix where
#' the columns correspond to independent replicates of observations.
#' @param A An observation matrix that links the measurement location to the finite elemen basis.
#' @param sigma.e The standard deviation of the measurement noise.
#'
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
#' Y = as.vector(A\%*\%u + sigma.e*rnorm(10))
#'
#' #compute log-likelihood of the data
#' lik1 <- rSPDE.loglike(op, Y, A, sigma.e)
#' cat(lik1)

rSPDE.loglike <- function(obj, Y, A, sigma.e)
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
  mu.post <- solve(R.post, AtY, system = "A")
  
  lik = n.rep * (prior.ld - posterior.ld - dim(A)[1] * log(2*pi) - sum(log(nugget))) / 2
  
  if (n.rep > 1) {
    lik = lik - 0.5 * sum(colSums(mu.post * (obj$Q %*% mu.post)))
    v = Q.e%*%(Y - A %*% mu.post)
    lik = lik - 0.5 * sum(colSums((Y - A %*% mu.post)*v)) 
  } else {
    lik = lik - 0.5 * (t(mu.post) %*% obj$Q %*% mu.post + t(Y - A %*% mu.post) %*% Q.e %*% (Y - A %*% mu.post))
  }
  return(as.double(lik))
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
#' Y = as.matrix(A\%*\%u + sigma.e*noise)
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
#' Y = as.matrix(A\%*\%u + sigma.e*noise)
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
