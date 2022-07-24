#'  Rational approximations of fractional SPDEs.
#'
#' \code{rSPDE} is used for approximating fractional elliptic SPDEs
#' \deqn{L^\beta (\tau u(s)) = W,}
#' where \eqn{L} is a differential operator and \eqn{\beta>0}
#' is a general fractional power.
#'
#' The approximation is based on a rational approximation of the
#' fractional operator, and allows for computationally efficient
#' inference and simulation.
#'
#' The main functions for computing rational approximation objects are:
#' \itemize{
#' \item{\code{\link{fractional.operators}}}{works for general
#' rational operators}
#' \item{\code{\link{matern.operators}}}{ works for random fields with
#' stationary Matern covariance functions}
#' \item{\code{\link{spde.matern.operators}}}{ works for random fields with
#' defined as solutions to a possibly non-stationary Matern-type SPDE model.}
#' \item{\code{\link{rspde.matern}}} {R-INLA implementation of the
#' covariance-based rational approximation for random fields with
#' stationary Matern covariance functions}
#' }
#' Basic statistical operations such as likelihood evaluations (see
#' \code{\link{rSPDE.loglike}, \link{rSPDE.matern.loglike}}) and kriging
#' predictions (see \code{\link{predict.rSPDEobj}, \link{predict.CBrSPDEobj}})
#' using the rational approximations are also implemented.
#'
#' For illustration purposes, the package contains a simple FEM implementation
#' for models on R. For spatial models,
#' the FEM implementation in the \code{R-INLA} package is recommended.
#'
#' For a more detailed introduction to the package, see the rSPDE Vignettes.
#'
#' @docType package
#' @name rSPDE
#' @import Matrix
#' @importFrom stats rnorm approx
NULL
