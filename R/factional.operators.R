

#' Rational approximations of fractional operators
#'
#' \code{fractional.operators} is used for computing an approximation, which can be used for inference and
#' simulation, of the fractional SPDE
#' \deqn{L^\beta (\tau u(s)) = W.}
#' Here \eqn{L} is a differential operator, \eqn{\beta>0} is
#' the fractional power, \eqn{\tau>0} is a scalar that determines the variance of the solution \eqn{u}, and
#' \eqn{W} is white noise.
#'
#' @param L A finite element discretization of the operator \eqn{cL}{cL}.
#' @param beta The positive fractional power.
#' @param C The mass matrix of the finite element discretization.
#' @param scale.factor A constant \eqn{c} so that the smallest eigenvalue of the non-discretized
#' operator \eqn{latex}{cL} is one.
#' @param m The order of the rational approximation, which needs to be a positive integer. The default value is 1.
#' Higer values gives a more accurate approximation, which are more computationally expensive to use for inference.
#' Currently, the largest value of m that is implemented is 4.
#' @param tau The constant that scales the variance of the solution. The default value is 1.
#'
#' @return \code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the
#' following arguments:
#' \item{Pl}{The operator \eqn{P_l}.}
#' \item{Pr}{The operator \eqn{P_r}.}
#' \item{C}{The mass matrix.}
#' \item{m}{The order of the rational approximation.}
#' \item{beta}{The fractional power.}
#' \item{type}{String indicating the type of approximation.}
#' \item{Q}{The matrix \code{t(Pl)\%*\%solve(C,Pl)}.}
#' \item{type}{String indicating the type of approximation.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso \code{\link{matern.operators}}, \code{\link{spde.matern.operators}}
#' @details The approximation is based on a rational approximation of the fractional operator,
#' resulting in an
#' approximate model on the form \deqn{P_l u(s) = P_r W,}
#' where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms of polynomials \eqn{p_j} for
#' \eqn{j=l,r}. The order of \eqn{p_r} is given by \code{m} and the order of \eqn{p_l} is \eqn{m + m_\beta}
#' where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and
#' \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})}
#' and \eqn{Q = P_l^T C^{-1} P_l}.
#'
#' @examples
#' #Compute rational approximation of a Gaussian process with Matern covariance function on R
#' kappa = 10
#' sigma = 1
#' nu = 0.8
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation of covariance function at 0.5
#' tau = sqrt(gamma(nu)/(sigma^2*kappa^(2*nu)*(4*pi)^(1/2)*gamma(nu+1/2)))
#' op <- fractional.operators(L = fem$G/kappa^2 + fem$C, beta=(nu + 1/2)/2,
#'                            C=fem$C, scale.factor = kappa^2, tau = tau)
#'
#' v = rep(0,101);v[51] = 1
#' c.approx = op$Pr \%*\% solve(op$Q, op$Pr \%*\% v)
#'
#' #plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x-0.5), kappa, nu, sigma), type = "l", ylab = "C(h)",
#'      xlab="h", main = "Matern covariance and rational approximations")
#' lines(x, c.approx, col = 2)

fractional.operators <- function(L,
                                 beta,
                                 C,
                                 scale.factor,
                                 m = 1,
                                 tau = 1)
  {
    C = Matrix::Diagonal(dim(C)[1], rowSums(C))
    Ci = Matrix::Diagonal(dim(C)[1], 1 / rowSums(C))
    I = Matrix::Diagonal(dim(C)[1])
    if (beta %% 1 == 0) {
      #not fractional
      Pr <- I
      Pl <- L
      if (beta > 1) {
        for (i in 1:beta) {
          Pl <- Pl %*% Ci %*% L
        }
      }

      Pl <- tau * Pl * scale.factor ^ beta
    } else {
      roots <- get.roots(m, beta)
      m_beta = max(1, floor(beta))

      CiL <- Ci %*% L

      #construct Pl
      Pl = I - CiL * roots$rb[1]
      if (length(roots$rb) > 1)
      {
        for (i in 2:length(roots$rb))
        {
          Pl <- Pl %*% (I - CiL * roots$rb[i])
        }
      }

      Lp <- C
      if (m_beta > 1) {
        for (i in 1:(m_beta - 1)) {
          Lp <- Lp %*% CiL
        }
      }
      Pl <- tau * Lp %*% Pl * scale.factor ^ beta / roots$factor

      #construct Pr
      Pr <- I - CiL * roots$rc[1]
      if (length(roots$rc) > 1) {
        for (i in 2:length(roots$rc)) {
          Pr <- Pr %*% (I - CiL * roots$rc[i])
        }
      }
    }

    Q = t(Pl) %*% Ci %*% Pl

    output <- list(Q = Q,
                   Pl = Pl,
                   Pr = Pr,
                   Ci = Ci,
                   C = C,
                   m = m,
                   beta = beta,
                   type = "fractional approximation")
    class(output) <- "rSPDEobj"
    return(output)
  }

#' Rational approximations of stationary Gaussian Matern random fields
#'
#' \code{matern.operators} is used for computing a rational SPDE approximation of a stationary Gaussian random
#' fields on \eqn{\mathbb{R}^d}{R^d} with a Matern covariance function
#' \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h)}{C(h) =
#' (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)}
#'
#' @param kappa Range parameter of the covariance function.
#' @param sigma Standard deviation of the covariance function.
#' @param nu Shape parameter of the covariance function.
#' @param G The stiffness matrix of a finite element discretization of the domain of interest.
#' @param C The mass matrix of a finite element discretization of the domain of interest.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a positive integer.
#' The default value is 1.
#'
#' @details The approximation is based on a rational approximation of the fractional operator
#' \eqn{(\kappa^2 -\Delta)^\beta}, where \eqn{\beta = (\nu + d/2)/2}.
#' This results in an approximate model on the form \deqn{P_l u(s) = P_r W,}
#' where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms of polynomials \eqn{p_j} for
#' \eqn{j=l,r}. The order of \eqn{p_r} is given by \code{m} and the order of \eqn{p_l} is \eqn{m + m_\beta}
#' where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and
#' \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})}
#' and \eqn{Q = P_l^T C^{-1} P_l}.
#'
#' @return \code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the
#' following arguments:
#' \item{Pl}{The operator \eqn{P_l}.}
#' \item{Pr}{The operator \eqn{P_r}.}
#' \item{C}{The mass matrix.}
#' \item{m}{The order of the rational approximation.}
#' \item{beta}{The fractional power.}
#' \item{Q}{The matrix \code{t(Pl)\%*\%solve(C,Pl)}.}
#' \item{type}{String indicating the type of approximation.}
#' \item{kappa}{Range parameter.}
#' \item{sigma}{Standard deviation.}
#' \item{nu}{Shape parameter.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso \code{\link{fractional.operators}}, \code{\link{spde.matern.operators}}
#'
#' @examples
#' #Compute rational approximation of a Gaussian process with Matern covariance function on R
#' kappa = 10
#' sigma = 1
#' nu = 0.8
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#' #compute rational approximation of covariance function at 0.5
#' op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
#'                        G = fem$G, C = fem$C, d = 1)
#' v = rep(0,101);v[51] = 1
#' c.approx = op$Pr \%*\% solve(op$Q, op$Pr \%*\% v)
#'
#' #plot the result and compare with the true Matern covariance
#' plot(x, matern.covariance(abs(x-0.5), kappa, nu, sigma), type = "l", ylab = "C(h)",
#'      xlab="h", main = "Matern covariance and rational approximation")
#' lines(x,c.approx,col=2)

matern.operators <- function(kappa,
                             sigma,
                             nu,
                             G,
                             C,
                             d = 2,
                             m = 1)
{
  tau = sqrt(gamma(nu) / (sigma^2 * kappa^(2*nu) * (4*pi)^(d /2) * gamma(nu + d/2)))
  beta = (nu + d/2)/2
  operators <- fractional.operators(L = G/kappa^2 + C,
                                    beta = beta,
                                    C = C,
                                    scale.factor = kappa^2,
                                    m = m,
                                    tau = tau)
  output <- operators
  output$kappa = kappa
  output$sigma = sigma
  output$nu = nu
  output$type = "Matern approximation"
  return(output)
}

#' Rational approximations of stationary Gaussian SPDE Matern random fields
#'
#' \code{spde.matern.operators} is used for computing a rational SPDE approximation of a Gaussian random
#' fields on \eqn{\mathbb{R}^d}{R^d} defined as a solution to the SPDE
#' \deqn{(\kappa(s) - \Delta)^\beta (\tau(s)u(s)) = W}
#'
#' @param kappa Vector with the, possibly spatially varying, range parameter evaluated at the locations
#' of the mesh used for the finite element discretization of the SPDE.
#' @param tau Vector with the, possibly spatially varying, precision parameter evaluated at the locations
#' of the mesh used for the finite element discretization of the SPDE.
#' @param nu Shape parameter of the covariance function, related to \eqn{\beta} through the equation
#' \eqn{\beta = (\nu + d/2)/2}.
#' @param G The stiffness matrix of a finite element discretization of the domain of interest.
#' @param C The mass matrix of a finite element discretization of the domain of interest.
#' @param d The dimension of the domain. The default value is 2.
#' @param m The order of the rational approximation, which needs to be a positive integer.
#' The default value is 1.
#'
#'
#' @details The approximation is based on a rational approximation of the fractional operator
#' \eqn{(\kappa(s)^2 -\Delta)^\beta}, where \eqn{\beta = (\nu + d/2)/2}.
#' This results in an approximate model on the form \deqn{P_l u(s) = P_r W,}
#' where \eqn{P_j = p_j(L)} are non-fractional operators defined in terms of polynomials \eqn{p_j} for
#' \eqn{j=l,r}. The order of \eqn{p_r} is given by \code{m} and the order of \eqn{p_l} is \eqn{m + m_\beta}
#' where \eqn{m_\beta} is the integer part of \eqn{\beta} if \eqn{\beta>1} and \eqn{m_\beta = 1} otherwise.
#'
#' The discrete approximation can be written as \eqn{u = P_r x} where \eqn{x \sim N(0,Q^{-1})}{x ~ N(0,Q^{-1})}
#' and \eqn{Q = P_l^T C^{-1} P_l}.
#'
#' @return \code{fractional.operators} returns an object of class "rSPDEobj". This is a list that contains the
#' following arguments:
#' \item{Pl}{The operator \eqn{P_l}.}
#' \item{Pr}{The operator \eqn{P_r}.}
#' \item{C}{The mass matrix.}
#' \item{m}{The order of the rational approximation.}
#' \item{beta}{The fractional power.}
#' \item{type}{String indicating the type of approximation.}
#' \item{Q}{The matrix \code{t(Pl)\%*\%solve(C,Pl)}.}
#' \item{nu}{Shape parameter.}
#' @export
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @seealso \code{\link{fractional.operators}}, \code{\link{spde.matern.operators}}
#'
#' @examples
#' #Sample non-stationary Matern field on R
#' tau = 1
#' nu = 0.8
#'
#' #create mass and stiffness matrices for a FEM discretization
#' x = seq(from = 0, to = 1, length.out = 101)
#' fem <- rSPDE.fem1d(x)
#'
#'  #define a non-stationary range parameter
#'  kappa <- seq(from = 2, to = 20, length.out = length(x))
#'
#' #compute rational approximation
#' op <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu,
#'                             G = fem$G, C = fem$C, d = 1)
#'
#' #sample the field
#' u <- simulate(op)
#'
#' #plot the sample
#' plot(x, u, type = "l", ylab = "u(s)", xlab = "s")

spde.matern.operators <- function(kappa,
                                  tau,
                                  nu,
                                  G,
                                  C,
                                  d = 2,
                                  m = 1)
{
  beta <- (nu + d/2)/2
  kp <- Matrix::Diagonal(dim(C)[1], kappa^2)
  Phi <- Matrix::Diagonal(dim(C)[1], 1/tau)
  scale.factor <- min(kappa)^2
  operators <- fractional.operators(L = (G + C %*% kp)/scale.factor,
                                    beta = beta,
                                    C = C,
                                    scale.factor = scale.factor,
                                    m = m,
                                    tau = 1)
  operators$Pr <- Phi %*% operators$Pr
  output <- operators
  output$beta = beta
  output$type = "Matern SPDE approximation"
  return(output)
}
