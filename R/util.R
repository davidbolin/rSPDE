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
#' \deqn{C(h) = \frac{\sigma^2}{2^(\nu-1)\Gamma(\nu)}(\kappa h)^\nu K_\nu(\kappa h)}{C(h) =
#' (\sigma^2/(2^(\nu-1)\Gamma(\nu))(\kappa h)^\nu K_\nu(\kappa h)}
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

matern.covariance <- function(h, kappa, nu, sigma)
{
  if (nu == 1/2) {
    C = sigma^2*exp(-kappa*abs(h))
  } else {
    C = (sigma^2 / (2^(nu - 1) * gamma(nu))) * ((kappa*abs(h))^nu) * besselK(kappa*abs(h), nu)
  }
  C[h == 0] = sigma^2
  return(C)
}

#' Summarise excurobj objects
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
  out$n = dim(object$L2)[1]
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
#' This function computes mass and stiffness matrices for a FEM approximation on R.
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
