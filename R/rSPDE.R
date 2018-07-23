#'  Rational approximations of fractional SPDEs.
#'
#' \code{rSPDE} is used for approximating fractional elliptic SPDEs \eqn{$L^\beta u(s) = W$}{L^\betau(s) = W},
#' where \eqn{L} is a differential operator and \eqn{\beta>0} is a general fractional power.
#'
#' The approximation is based on a rational approximation of the fractional operator, and allows for computationally
#' efficient inference and simulation.
#'
#' The main function for computing the rational operators is \code{\link{fractional.operators}}, and a simplified
#' interface for the popular case of Matern covariances is provided through the function
#' \code{\link{matern.operators}}. Basic statistical operations such as likelihood evaluations
#' (see \code{\link{rSPDE.loglike}}) and kriging predictions (see \code{\link{rSPDE.krig}}) using the fractional
#' approximations are also implemented.
#'
#' For illustration purposes, the package contains a simple FEM implementation for models on R. For spatial models,
#' the FEM implementation in the \code{R-INLA} package is recommended.
#'
#' @docType package
#' @name rSPDE
#' @import Matrix
NULL
