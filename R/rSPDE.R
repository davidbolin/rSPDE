#'  Rational approximations of fractional SPDEs.
#'
#' \code{rSPDE} is used for approximating fractional elliptic SPDEs
#' \eqn{$L^\beta u(s) = W$}{L^\betau(s) = W},
#' where \eqn{L} is a differential operator and \eqn{\beta>0} is a general fractional power.
#'
#' The approximation is based on a rational approximation of the fractional operator,
#' and allows for computationally efficient inference and simulation.
#'
#' The main function for computing the rational operators is \code{\link{fractional.operators}},
#' and the following simplified interfaces are available
#' \itemize{
#' \item{\code{\link{matern.operators}}}{Compuation of operators for random fields with
#' stationary Matern covariance functions}
#' \item{\code{\link{spde.matern.operators}}}{Compuation of operators for random fields with
#' defined as solutions to a possibly non-stationary Matern-type SPDE model.}
#' }
#' Basic statistical operations such as likelihood evaluations (see \code{\link{rSPDE.loglike}}) and kriging
#' predictions (see \code{\link{predict.rSPDEobj}}) using the fractional approximations are also implemented.
#'
#' For illustration purposes, the package contains a simple FEM implementation for models on R. For spatial models,
#' the FEM implementation in the \code{R-INLA} package is recommended.
#' 
#' For a more detailed introduction to the package, see the rSPDE Vignette.
#'
#' @docType package
#' @name rSPDE
#' @import Matrix
#' @importFrom stats rnorm approx
NULL
