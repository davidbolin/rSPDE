#'  Rational approximations of fractional SPDEs.
#'
#' `rSPDE` is used for approximating fractional elliptic SPDEs
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
#' \item{[fractional.operators()]}{works for general
#' rational operators}
#' \item{[matern.operators()]}{ works for random fields with
#' stationary Matern covariance functions}
#' \item{[spde.matern.operators()]}{ works for random fields with
#' defined as solutions to a possibly non-stationary Matern-type SPDE model.}
#' \item{[rspde.matern()]} {R-INLA implementation of the
#' covariance-based rational approximation for random fields with
#' stationary Matern covariance functions}
#' }
#' Basic statistical operations such as likelihood evaluations (see
#' `[rSPDE.loglike], [rSPDE.matern.loglike]`) and kriging
#' predictions (see `[predict.rSPDEobj], [predict.CBrSPDEobj]`)
#' using the rational approximations are also implemented.
#'
#' For illustration purposes, the package contains a simple FEM implementation
#' for models on R. For spatial models,
#' the FEM implementation in the `R-INLA` package is recommended.
#'
#' For a more detailed introduction to the package, see the rSPDE Vignettes.
#'
#' @docType package
#' @name rSPDE
#' @aliases rSPDE-package
#' @import Matrix
#' @importFrom stats rnorm approx quantile
#' @importFrom stats dnorm pnorm dbeta nobs deviance logLik
#' @importFrom methods as
#' @importFrom stats simulate lm logLik na.omit optim sd terms
#' @importFrom fmesher fm_mesh_2d fm_basis fm_block fm_row_kron fm_mesh_1d fm_fem
#' @importFrom lifecycle deprecated
#' @importFrom broom augment glance
#' @export augment glance
#'
NULL

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("rSPDE")
   packageStartupMessage(
    "This is rSPDE ", version, "\n",
    "- See https://davidbolin.github.io/rSPDE for vignettes and manuals."
  )
}