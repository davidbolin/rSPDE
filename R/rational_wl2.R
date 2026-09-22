## Weighted-L2 rational coefficients with spectral interval.
##
## Rational approximations of x^target on [x_min, 1] (x = 1/lambda, lambda the
## eigenvalues of L_h = kappa^{-2}(kappa^2 C + G) relative to C, so lambda >= 1)
## that minimise the weighted L2 error with the Weyl weight
##     w_d(x) = x^(-1-d/2) (1-x)^(d/2-1).
##
## The fit is done by variable projection: the poles are the outer parameters
## (optimised by Levenberg-Marquardt) and the residues are obtained by
## non-negative least squares, so that every returned set of coefficients
## defines a valid (positive definite) model.
#

#' @name wl2_weyl_cells
#' @title Nodes and cell integrals of the Weyl weight
#' @description Log-equispaced nodes on \[x_min, 1\] together with the integrals
#' of the Weyl density over the corresponding cells, computed with a
#' sub-midpoint rule in u = log(x). The weight is never evaluated at x = 1,
#' where it is singular for d = 1.
#' @param d Dimension.
#' @param x_min Lower end of the interval.
#' @param n_grid Number of nodes.
#' @param n_sub Number of sub-intervals used for each cell integral.
#' @return A list with elements `x` (nodes) and `w` (cell integrals).
#' @noRd
wl2_weyl_cells <- function(d, x_min, n_grid = 300, n_sub = 40, s = 0) {
  lx <- seq(log(x_min), 0, length.out = n_grid)
  ed <- c(lx[1], (lx[-1] + lx[-n_grid]) / 2, 0)
  w <- numeric(n_grid)
  for (i in seq_len(n_grid)) {
    e <- seq(ed[i], ed[i + 1], length.out = n_sub + 1)
    xm <- exp((e[-1] + e[-(n_sub + 1)]) / 2)
    ## The integral is taken in log x, so the Weyl density w_d(x) appears
    ## multiplied by x; the data weight x^s shifts the same exponent.
    w[i] <- sum(xm^(s - d / 2) * (1 - xm)^(d / 2 - 1)) * (e[2] - e[1])
  }
  list(x = exp(lx), w = w)
}

#' @name wl2_ls
#' @title Least squares solution tolerating rank deficiency
#' @description Aliased columns get coefficient zero instead of `NA`. The rank
#' tolerance has to be far tighter than the default of `.lm.fit()`: the columns
#' span poles over six orders of magnitude and are strongly correlated, and at
#' the default 1e-7 the last columns are declared aliased and their residues
#' set to zero, which silently throws away terms of the approximation.
#'
#' @param A Design matrix.
#' @param b Right hand side.
#' @return The least squares coefficients.
#' @noRd
wl2_ls <- function(A, b) {
  cf <- .lm.fit(A, b, tol = 1e-14)$coefficients
  cf[!is.finite(cf)] <- 0
  cf
}

#' @name wl2_nnls
#' @title Non-negative least squares 
#' @description Minimises `||A x - b||` subject to `x >= 0`. Written for the
#' small problems (at most eight columns) arising in the variable projection,
#' so that no external NNLS dependency is needed.
#' @param A Design matrix.
#' @param b Right hand side.
#' @return The non-negative least squares solution.
#' @noRd
wl2_nnls <- function(A, b) {
  n <- ncol(A)
  if (n == 0) {
    return(numeric(0))
  }
  ## The columns span several orders of magnitude, because the poles do: for
  ## m = 6 the norm of the last column is a millionth of the first. The
  ## optimality test compares gradients across columns, so it is only
  ## meaningful once the columns are on a common scale. Without this the
  ## algorithm stops before the last columns are brought in, and silently
  ## drops terms of the approximation. The scaling is positive and diagonal,
  ## so the non-negativity constraints are unchanged.
  scl <- sqrt(base::colSums(A^2))
  scl[scl <= 0] <- 1
  A <- A / rep(scl, each = nrow(A))
  x <- numeric(n)
  P <- logical(n)
  w <- drop(base::crossprod(A, b))
  tol <- 10 * .Machine$double.eps * max(dim(A)) * max(abs(w))
  it <- 0
  maxit <- 10 * n
  while (any(!P) && max(w[!P]) > tol && it < maxit) {
    it <- it + 1
    idx <- which(!P)
    P[idx[which.max(w[idx])]] <- TRUE
    inner <- 0
    repeat {
      inner <- inner + 1
      s <- numeric(n)
      s[P] <- wl2_ls(A[, P, drop = FALSE], b)
      if (min(s[P]) >= 0) {
        x <- s
        break
      }
      neg <- P & (s < 0)
      ratios <- x[neg] / (x[neg] - s[neg])
      step <- min(ratios[is.finite(ratios)], 1)
      x <- x + step * (s - x)
      P <- P & (x > 10 * .Machine$double.eps * max(abs(x)))
      x[!P] <- 0
      if (!any(P) || inner > maxit) {
        break
      }
    }
    w <- drop(base::crossprod(A, b - base::`%*%`(A, x)))
  }
  x / scl
}

#' @name wl2_resjac
#' @title Residual and Golub-Pereyra Jacobian of the variable projection
#' @description The residual is `r(theta) = A(theta) c - f` with `c` the
#' non-negative least squares solution. On the passive set P the coefficients
#' are `c_P = A_P^+ f`, so `r(theta) = -(I - P_A) f` and
#' \deqn{\partial r/\partial\theta_k = (I - P_A) (\partial A_P/\partial\theta_k) c_P
#'   + (A_P^+)^T (\partial A_P/\partial\theta_k)^T (-r),}
#' the Golub-Pereyra derivative. With `u_j = x/(1 + (b_j - 1)x)` and
#' `b_j = exp(theta_j)` the column derivative is `-b_j u_j^2`, and for the
#' shifted classes the integer factor `g^q` with
#' `g = x/(1 + (b_0 - 1)x)` and `theta_1 = logit(b_0)` moves every column
#' through `d(g^q)/dtheta_1 = -q g^(q+1) b_0 (1 - b_0)`. Columns whose residue
#' is zero do not enter, so their derivative is zero.
#' @param kind `"plain"`, `"shifted"` or `"shifted2"`; see `wl2_q()`.
#' @param x,sw Quadrature nodes and the square roots of the cell integrals.
#' @param f The weighted target.
#' @param theta The outer parameters.
#' @param need_jac Compute the Jacobian, or only the residual?
#' @return A list with `r`, `c` and, if requested, `J`.
#' @noRd
wl2_resjac <- function(kind, x, sw, f, theta, need_jac = TRUE,
                       k_term = FALSE) {
  mm <- base::`%*%`
  q <- wl2_q(kind)
  shifted <- q > 0
  b <- exp(pmin(pmax(if (shifted) theta[-1] else theta, -5), 60))
  n <- length(x)
  M <- length(b)
  U <- matrix(x / (1 + rep(b - 1, each = n) * x), n, M)
  if (shifted) {
    b0 <- 1 / (1 + exp(-theta[1]))
    g <- x / (1 + (b0 - 1) * x)
    A <- (g^q * U) * sw
  } else {
    A <- U * sw
  }
  if (k_term) {
    A <- cbind(A, sw)
  }
  cc <- wl2_nnls(A, f)
  r <- drop(mm(A, cc)) - f
  if (!need_jac) {
    return(list(r = r, c = cc))
  }
  keep <- cc > 0
  if (!any(keep)) {
    return(list(r = r, c = cc, J = matrix(0, n, length(theta))))
  }
  ## base::qr, not the S4 generic of Matrix: these are plain matrices and the
  ## dispatch costs more than the factorisation. Q is formed once and reused
  ## for every parameter; going through qr.resid()/qr.qy() per parameter
  ## instead is 2.7 times slower, since each does a full pass over the
  ## factorisation.
  qrA <- base::qr(A[, keep, drop = FALSE])
  Q <- base::qr.Q(qrA)
  R <- base::qr.R(qrA)
  J <- matrix(0, n, length(theta))
  for (k in seq_along(theta)) {
    D <- matrix(0, n, ncol(A))
    if (shifted && k == 1) {
      ## d(g^q)/dtheta_0 = q g^(q-1) dg/dtheta_0 = -q g^(q+1) b_0 (1 - b_0)
      D[, seq_len(M)] <- ((-q * g^(q + 1) * b0 * (1 - b0)) * U) * sw
    } else {
      j <- if (shifted) k - 1 else k
      D[, j] <- (if (shifted) g^q else 1) * (-b[j] * U[, j]^2) * sw
    }
    v <- drop(mm(D, cc))
    t1 <- v - mm(Q, base::crossprod(Q, v))
    z <- drop(base::crossprod(D[, keep, drop = FALSE], -r))
    J[, k] <- t1 + mm(Q, backsolve(R, z, transpose = TRUE))
  }
  list(r = r, c = cc, J = J)
}

#' @name wl2_lm
#' @title Levenberg-Marquardt minimisation of a residual vector
#' @description Minimises `sum(res(theta)^2)` with a forward-difference
#' Jacobian. Used for the outer (pole) problem, which has at most eight
#' parameters.
#' @param res Function returning the residual vector.
#' @param theta Starting value.
#' @param ftol,xtol Relative tolerances on the objective and on the step.
#' @param maxfev Maximum number of evaluations of `res`.
#' @param maxiter Maximum number of Levenberg-Marquardt iterations.
#' @return A list with elements `theta`, `cost` (half the sum of squares) and
#' `nfev`.
#' @noRd
wl2_lm <- function(res, theta, ftol = wl2_tol, xtol = wl2_tol,
                   maxfev = 3000, maxiter = 400) {
  np <- length(theta)
  ev <- res(theta, TRUE)
  analytic <- !is.null(ev$J)
  r <- ev$r
  f <- sum(r^2)
  nfev <- 1
  lambda <- 1e-3
  J <- matrix(0, length(r), np)
  for (iter in seq_len(maxiter)) {
    if (analytic) {
      J <- ev$J
    } else {
      h <- sqrt(.Machine$double.eps) * pmax(abs(theta), 1)
      for (j in seq_len(np)) {
        tj <- theta
        tj[j] <- tj[j] + h[j]
        J[, j] <- (res(tj, FALSE)$r - r) / h[j]
        nfev <- nfev + 1
      }
    }
    if (nfev > maxfev) {
      break
    }
    JtJ <- crossprod(J)
    Jtr <- as.vector(crossprod(J, r))
    dg <- diag(JtJ)
    dg <- pmax(dg, 1e-10 * max(dg, 1e-300))
    accepted <- FALSE
    while (lambda < 1e14) {
      step <- tryCatch(
        as.vector(solve(JtJ + lambda * diag(dg, np), -Jtr)),
        error = function(e) NULL
      )
      if (is.null(step) || any(!is.finite(step))) {
        lambda <- lambda * 10
        next
      }
      theta_new <- theta + step
      ev_new <- res(theta_new, analytic)
      r_new <- ev_new$r
      f_new <- sum(r_new^2)
      nfev <- nfev + 1
      if (is.finite(f_new) && f_new < f) {
        df <- f - f_new
        dx <- sqrt(sum(step^2))
        theta <- theta_new
        ev <- ev_new
        r <- r_new
        f <- f_new
        lambda <- max(lambda / 10, 1e-12)
        accepted <- TRUE
        if (df <= ftol * f || dx <= xtol * (sqrt(sum(theta^2)) + xtol)) {
          return(list(theta = theta, cost = f / 2, nfev = nfev))
        }
        break
      }
      lambda <- lambda * 10
      if (nfev > maxfev) {
        break
      }
    }
    if (!accepted || nfev > maxfev) {
      break
    }
  }
  list(theta = theta, cost = f / 2, nfev = nfev)
}

#' @name wl2_tol
#' @title Convergence tolerance of the outer minimisation
#' @description Relative tolerance on the objective and on the step in the
#' Levenberg-Marquardt loop. The objective is flat near its minimum, so the
#' last digits of the poles are worth nothing in the error they produce:
#' loosening this from 1e-13 to 1e-10 leaves the weighted-L2 error of every
#' fit in the tables unchanged to five significant digits, and the reference
#' coefficients unchanged to six, while cutting the cost of building a table
#' by a third.
#' @noRd
wl2_tol <- 1e-10

#' @name wl2_have_cpp
#' @title Is the compiled weighted-L2 solver available?
#' @description The Levenberg-Marquardt loop, the variable projection and the
#' non-negative least squares are also implemented in C (`src/wl2_fit.cpp`),
#' which is part of every build. The R implementation is kept as a reference
#' and as a fallback, and is selected by setting the option
#' `rSPDE.wl2.use.cpp` to `FALSE`.
#' @return `TRUE` if the compiled routine should be used.
#' @noRd
wl2_have_cpp <- function() {
  if (!isTRUE(getOption("rSPDE.wl2.use.cpp", TRUE))) {
    return(FALSE)
  }
  if (is.null(wl2_cache$have_cpp)) {
    wl2_cache$have_cpp <- is.character(tryCatch(
      getNativeSymbolInfo("rspde_wl2_lm", PACKAGE = "rSPDE")$name,
      error = function(e) NULL
    ))
  }
  wl2_cache$have_cpp
}

#' @name wl2_lm_run
#' @title Run the outer minimisation, in C if possible
#' @description Dispatches to the compiled solver, falling back to the R
#' implementation if it is not available or if it fails. The two give the same
#' answer to roundoff; see the tests.
#' @param kind One of "plain", "shifted" or "shifted2"; see `wl2_q()`.
#' @param x,sw,f Nodes, square roots of the cell integrals, weighted target.
#' @param theta0 Starting value.
#' @param tol Relative tolerance on the objective and on the step; see
#' `wl2_tol`.
#' @return A list with `theta` and `cost`, as `wl2_lm()`.
#' @noRd
wl2_lm_run <- function(kind, x, sw, f, theta0, tol = wl2_tol,
                       k_term = FALSE) {
  if (wl2_have_cpp()) {
    o <- tryCatch(
      .Call("rspde_wl2_lm", as.double(x), as.double(sw), as.double(f),
        as.double(theta0), as.integer(wl2_q(kind)), k_term,
        as.double(c(tol, tol, 3000, 400)),
        PACKAGE = "rSPDE"
      ),
      error = function(e) NULL
    )
    if (!is.null(o) && all(is.finite(o$theta)) && is.finite(o$cost)) {
      return(o)
    }
  }
  wl2_lm(function(theta, need_jac = TRUE) {
    wl2_resjac(kind, x, sw, f, theta, need_jac, k_term)
  }, theta0, ftol = tol, xtol = tol)
}

#' @name wl2_q
#' @title Power of the integer factor of a class
#' @description The covariance classes are
#' \eqn{r(x) = (x/(1-p_0x))^q \sum_j r_j x/(1-p_jx)} with `q` equal to
#' \eqn{\lfloor\alpha\rfloor}. `q = 0` is the plain class, which has no
#' shift parameter at all; `q = 1` and `q = 2` share one shift `p_0` across the
#' whole factor. A separate shift per factor was tried and is not worth it: at
#' `q = 2` two free shifts match one shared shift to between 0.98 and 1.00 of
#' the error over the cases tested, the shared one being a stationary point of
#' the larger problem.
#' @param kind `"plain"`, `"shifted"` or `"shifted2"`.
#' @return 0, 1 or 2.
#' @noRd
wl2_q <- function(kind) {
  switch(kind, plain = 0L, shifted = 1L, shifted2 = 2L,
    stop("unknown class: ", kind)
  )
}

#' @name wl2_block_poly
#' @title Polynomial of a weighted-L2 precision block
#' @description The block of the precision matrix belonging to the pole
#' \eqn{p_i} is
#' \eqn{(L - p_0C)C^{-1}\cdots(L - p_0C)C^{-1}(L - p_iC)/r_i}, with the
#' shifted factor repeated `q` times. Writing \eqn{T = C^{-1}L} this is
#' \eqn{C(T-p_0)^q(T-p_i)/r_i}, so the block is a linear combination of the
#' matrices \eqn{P_j = CT^j}, `j` from 0 to `q + 1`. This returns the
#' coefficients of that combination, before the division by \eqn{r_i}.
#' @param q The power of the integer factor; see `wl2_q()`.
#' @param p0 The shared shift, unused when `q` is 0.
#' @param p_i The pole of the block.
#' @return A numeric vector of length `q + 2`, the coefficient of \eqn{P_j}
#' in position `j + 1`.
#' @noRd
wl2_block_poly <- function(q, p0, p_i) {
  cf <- 1
  ## Multiply by (T - p_0) q times, then once by (T - p_i). Coefficients are
  ## in increasing order of the power, so multiplying by (T - a) shifts and
  ## subtracts.
  for (l in seq_len(q)) {
    cf <- c(0, cf) - p0 * c(cf, 0)
  }
  c(0, cf) - p_i * c(cf, 0)
}

#' @name wl2_mass_powers
#' @title The matrices \eqn{CT^j} of a finite element discretisation
#' @description With \eqn{L = C + G/\kappa^2} and \eqn{T = C^{-1}L},
#' \eqn{P_j = CT^j = \sum_{l=0}^j \binom{j}{l} G_l/\kappa^{2l}}, where
#' \eqn{G_0 = C}, \eqn{G_1 = G} and \eqn{G_l = GC^{-1}G_{l-1}}. These are
#' exactly the matrices the finite element assembly provides, so the blocks
#' need no inverse of the mass matrix.
#' @param q The power of the integer factor; powers up to `q + 1` are returned.
#' @param kappa The range parameter.
#' @param Gl A function of one argument returning \eqn{G_l}; `l = 0` must give
#' the mass matrix.
#' @return A list of length `q + 2`, with \eqn{P_j} in position `j + 1`.
#' @noRd
wl2_mass_powers <- function(q, kappa, Gl) {
  lapply(0:(q + 1L), function(j) {
    out <- Gl(0L)
    for (l in seq_len(j)) {
      out <- out + choose(j, l) * Gl(l) / kappa^(2 * l)
    }
    out
  })
}

#' @name wl2_operator_powers
#' @title The matrices \eqn{CT^j} from the operator matrices themselves
#' @description As `wl2_mass_powers()`, but forming
#' \eqn{P_j = L(C^{-1}L)^{j-1}} by explicit products, for the paths where the
#' higher finite element matrices \eqn{G_2, G_3, \ldots} are not assembled.
#' @param q The power of the integer factor; powers up to `q + 1` are returned.
#' @param L The shifted stiffness matrix \eqn{C + G/\kappa^2}.
#' @param C0 The mass matrix.
#' @param Ci The inverse of the mass matrix, usually the lumped one.
#' @return A list of length `q + 2`, with \eqn{P_j} in position `j + 1`.
#' @noRd
wl2_operator_powers <- function(q, L, C0, Ci) {
  P <- vector("list", q + 2L)
  P[[1]] <- C0
  P[[2]] <- L
  for (j in seq_len(q)) {
    P[[j + 2L]] <- P[[j + 1L]] %*% Ci %*% L
  }
  P
}

#' @name wl2_block_builder
#' @title Precision blocks of a weighted-L2 approximation
#' @description Combines the powers \eqn{P_j} of `wl2_mass_powers()` or
#' `wl2_operator_powers()` with the coefficients of `wl2_block_poly()` into a
#' function of the block index.
#' @param q The power of the integer factor; see `wl2_q()`.
#' @param P The powers, as returned by `wl2_mass_powers()`.
#' @param r,p The residues and poles.
#' @param p0 The shared shift, `NULL` when `q` is 0.
#' @return A function of `i` returning the `i`th block.
#' @noRd
wl2_block_builder <- function(q, P, r, p, p0) {
  function(i) {
    cf <- wl2_block_poly(q, p0, p[i])
    out <- cf[1] * P[[1]]
    for (j in seq_along(cf)[-1]) {
      out <- out + cf[j] * P[[j]]
    }
    out / r[i]
  }
}

#' @name wl2_basis
#' @title Basis matrix of the rational class
#' @param kind `"plain"`, `"shifted"` or `"shifted2"`; see `wl2_q()`.
#' @param x Nodes.
#' @param theta Outer parameters: log(b) for "plain", and logit(b_0) followed by
#' log(b) for the shifted classes.
#' @return The basis matrix, with one column per term.
#' @noRd
wl2_basis <- function(kind, x, theta, k_term = FALSE) {
  q <- wl2_q(kind)
  b <- exp(pmin(pmax(if (q > 0) theta[-1] else theta, -5), 60))
  n <- length(x)
  A <- matrix(x / (1 + rep(b - 1, each = n) * x), n, length(b))
  if (q > 0) {
    b0 <- 1 / (1 + exp(-theta[1]))
    A <- (x / (1 + (b0 - 1) * x))^q * A
  }
  if (k_term) {
    ## The constant term is the last column, so that the residues keep their
    ## positions; it is non-negative like the rest.
    A <- cbind(A, rep(1, n))
  }
  A
}

#' @name wl2_starts
#' @title Starting values for the outer problem
#' @noRd
wl2_starts <- function(M, x_min, n_starts, kind, start, seed) {
  Lp <- min(log(1 / x_min), 2.2 * M + 6)
  inits <- list()
  if (!is.null(start)) {
    inits[[length(inits) + 1]] <- as.numeric(start)
  }
  inits[[length(inits) + 1]] <- seq(0.05 * Lp, Lp, length.out = M)
  if (n_starts > 1) {
    ## A local RNG stream, so that the coefficients are reproducible and the
    ## user's random number stream is left untouched.
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1)
    }
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    set.seed(seed)
    for (i in seq_len(n_starts - 1)) {
      inits[[length(inits) + 1]] <- sort(stats::runif(M, 0, Lp + 1))
    }
  }
  if (wl2_q(kind) > 0) {
    inits <- lapply(inits, function(t) if (length(t) == M) c(1.5, t) else t)
  }
  inits
}

#' Weighted-L2 rational coefficients
#'
#' Computes the coefficients of a rational approximation of \eqn{x^{\alpha}}
#' (covariance type) or \eqn{x^{\alpha/2}} (operator type) on the interval
#' \eqn{[x_{\min}, 1]}, minimising the \eqn{L_2} error weighted by the Weyl
#' density \eqn{w_d(x) = x^{-1-d/2}(1-x)^{d/2-1}}.
#'
#' Here \eqn{x = 1/\lambda}, where \eqn{\lambda} are the eigenvalues of
#' \eqn{L_h = \kappa^{-2}(\kappa^2 C + G)} relative to \eqn{C}, so that
#' \eqn{\lambda \ge 1} and \eqn{x \in (0, 1]}. The value \eqn{x_{\min}}
#' corresponds to the largest eigenvalue, i.e. to the resolution of the mesh;
#' `x_min = NULL` gives the mesh-free fit on \eqn{(0,1]}.
#'
#' The classes are parameterised so that the resulting model is guaranteed to
#' be valid (all residues are non-negative and all poles are smaller than one).
#' For the covariance type the class carries an integer factor of order
#' \eqn{q = \lfloor\alpha\rfloor}, with one shift \eqn{p_0} shared by its
#' factors:
#' \describe{
#' \item{covariance}{\eqn{r(x) = (x/(1-p_0 x))^q\sum_{j=1}^m r_j x/(1-p_j x)},
#' for \eqn{q = 0}, 1 and 2, the factor being absent when \eqn{q = 0}.}
#' \item{operator}{\eqn{r(x) = \sum_{i=1}^{m+1} r_i x/(1-p_i x)}, approximating
#' \eqn{x^{\alpha/2}}.}
#' }
#' In contrast to the tabulated `"brasil"`, `"chebfun"` and `"chebfunLB"`
#' coefficients, there is no constant term, so that the covariance-based models
#' built from these coefficients have \eqn{m} instead of \eqn{m+1} blocks.
#'
#' The fit is a variable projection: the poles are found by
#' Levenberg-Marquardt, with the residues eliminated by non-negative least
#' squares at every step. This inner loop is compiled (`src/wl2_fit.cpp`) and
#' is a few times faster than the equivalent R code, which is kept as a
#' reference implementation and is used instead when the option
#' `rSPDE.wl2.use.cpp` is set to `FALSE`.
#'
#' @param alpha The exponent, \eqn{\alpha = \nu + d/2}. For the covariance type,
#' \eqn{\lfloor \alpha\rfloor} must be 0, 1 or 2; for the operator type,
#' \eqn{\alpha} must be smaller than 2.
#' @param d The dimension of the domain.
#' @param m The order of the rational approximation.
#' @param x_min Lower end of the spectral interval, or `NULL` for the mesh-free
#' fit on \eqn{(0,1]}. See [rspde.xmin()].
#' @param type Either `"covariance"` or `"operator"`.
#' @param n_starts Number of starting values for the outer optimisation. Eight
#' starts are needed for `m = 4`; for `m` at most 3 a single start suffices.
#' @param n_grid Number of quadrature nodes.
#' @param start Optional starting value for the outer parameters, for instance
#' the `theta` element returned by a fit for a neighbouring value of `alpha`.
#' @param seed Seed of the local random number stream used for the random
#' starting values. The user's random number stream is not affected.
#' @param s Exponent of an extra factor \eqn{x^s} in the weight, on top of the
#' Weyl density. `s = 0`, the default, is the weight for which the objective is
#' the \eqn{L_2} error of the covariance. Larger `s` de-emphasises the small
#' eigenvalues, and corresponds to measuring the error of the solution operator
#' against data in \eqn{H^s} rather than in \eqn{L_2}.
#' @param k_term Should a constant term be included, as in the tabulated
#' coefficients? It is fitted, non-negative, and returned as `k`. The default
#' `FALSE` is what the covariance-based models need, since a constant term
#' makes the trace infinite. It only makes sense together with a finite
#' `x_min`: on \eqn{(0,1]} with `s` below \eqn{d/2} the objective is infinite
#' for any non-zero constant, and the fit drives it to zero.
#' @param continuation Should the fit use continuation in `m`, solving for
#' `1, ..., m` and starting each order from the poles of the previous one with
#' one pole inserted? This is needed beyond `m = 4`, where a plain multistart
#' tends to settle at the optimum of the next smaller order. Ignored when
#' `start` is supplied.
#'
#' @return A list with elements
#' \item{r}{The residues, all non-negative.}
#' \item{p}{The poles, all smaller than one.}
#' \item{p0}{The shift of the integer factor, or `NULL` if
#' \eqn{\lfloor\alpha\rfloor = 0}.}
#' \item{q}{The power of the integer factor, \eqn{\lfloor\alpha\rfloor} for
#' the covariance type and 0 for the operator type.}
#' \item{k}{The constant term: 0 unless `k_term` is `TRUE`.}
#' \item{rel_err}{The relative weighted \eqn{L_2} error of the fit.}
#' \item{x_min}{The lower end of the interval that was used.}
#' \item{theta}{The internal parameters, for warm starts.}
#' \item{kind}{`"plain"`, `"shifted"` or `"shifted2"`, for `q` equal to 0, 1
#' and 2.}
#' @export
#' @seealso [rspde.xmin()], [matern.operators()]
#' @examples
#' cf <- rational.coefficients.wl2(alpha = 0.75, d = 1, m = 2)
#' cf$rel_err
rational.coefficients.wl2 <- function(alpha, d, m, x_min = NULL,
                                      type = c("covariance", "operator"),
                                      n_starts = 8, n_grid = 300,
                                      start = NULL, seed = 1L,
                                      continuation = TRUE,
                                      s = 0, k_term = FALSE) {
  type <- match.arg(type)
  if (length(s) != 1 || !is.finite(s)) {
    stop("s must be a single finite number.")
  }
  k_term <- isTRUE(k_term)
  if (continuation && is.null(start) && m > 1) {
    ## Continuation in m. A plain multistart is not reliable beyond m = 4: the
    ## optimiser then tends to settle at the optimum of the next smaller order,
    ## with one residue exactly zero, which wastes a block. Starting from the
    ## poles of the (m-1)-fit with one pole inserted at each gap finds the
    ## global optimum in those cases, and reproduces the multistart result
    ## where that already found it.
    prev <- NULL
    for (mm in seq_len(m)) {
      best <- wl2_fit(
        alpha, d, mm, x_min, type, n_starts, n_grid, NULL, seed, s, k_term
      )
      if (!is.null(prev)) {
        for (st in wl2_insert_pole(prev$theta, !is.null(prev$p0))) {
          cf <- wl2_fit(
            alpha, d, mm, x_min, type, 1L, n_grid, st, seed, s, k_term
          )
          if (cf$rel_err < best$rel_err) {
            best <- cf
          }
        }
      }
      prev <- best
    }
    return(prev)
  }
  wl2_fit(alpha, d, m, x_min, type, n_starts, n_grid, start, seed, s, k_term)
}

#' @name wl2_insert_pole
#' @title Starting values obtained by adding one pole to a fit
#' @param theta The internal parameters of a fit of order m - 1.
#' @param shifted Does the class have a shift, that is, is `q` positive?
#' @return A list of starting values of order m, one per insertion point.
#' @noRd
wl2_insert_pole <- function(theta, shifted) {
  lb <- sort(if (shifted) theta[-1] else theta)
  n <- length(lb)
  gaps <- c(lb[1] - 2, if (n > 1) (lb[-1] + lb[-n]) / 2, lb[n] + 2)
  lapply(gaps, function(g) {
    nb <- sort(c(lb, g))
    if (shifted) c(theta[1], nb) else nb
  })
}

#' @name wl2_fit
#' @title A single weighted-L2 fit
#' @description The fit itself; [rational.coefficients.wl2()] adds the
#' continuation in m.
#' @noRd
wl2_fit <- function(alpha, d, m, x_min, type, n_starts, n_grid, start, seed,
                    s = 0, k_term = FALSE) {
  if (m < 1) {
    stop("m must be a positive integer!")
  }
  if (type == "operator") {
    target <- alpha / 2
    M <- m + 1
    kind <- "plain"
    if (target >= 1) {
      stop("For type = 'operator', alpha must be smaller than 2.")
    }
  } else {
    m_alpha <- floor(alpha)
    if (!(m_alpha %in% c(0, 1, 2))) {
      stop(paste0(
        "The weighted-L2 coefficients are only implemented for ",
        "floor(alpha) equal to 0, 1 or 2, but floor(alpha) = ", m_alpha, "."
      ))
    }
    target <- alpha
    M <- m
    kind <- c("plain", "shifted", "shifted2")[m_alpha + 1]
  }
  ## Near x = 0 the integrand is e(x)^2 x^(s - 1 - d/2) with e(x) ~ x^target,
  ## so the objective is finite when 2 * target + s > d / 2. With the data
  ## weight x^s this is weaker than the finite-trace condition s = 0.
  if (2 * target + s <= d / 2) {
    stop(paste0(
      "The fit requires 2 * target + s > d / 2, i.e. a finite objective; ",
      "here 2 * ", signif(target, 4), " + ", signif(s, 4), " <= ", d / 2, "."
    ))
  }

  if (is.null(x_min)) {
    x_min <- 1e-26
    n_grid <- max(n_grid, 900)
  }

  cells <- wl2_weyl_cells(d, x_min, n_grid, s = s)
  x <- cells$x
  sw <- sqrt(cells$w)
  f <- x^target * sw

  best <- NULL
  for (theta0 in wl2_starts(M, x_min, n_starts, kind, start, seed)) {
    o <- tryCatch(wl2_lm_run(kind, x, sw, f, theta0, k_term = k_term),
      error = function(e) NULL
    )
    if (!is.null(o) && (is.null(best) || o$cost < best$cost)) {
      best <- o
    }
  }
  if (is.null(best)) {
    stop("The rational approximation could not be computed.")
  }

  theta <- best$theta
  A <- wl2_basis(kind, x, theta, k_term) * sw
  cc <- wl2_nnls(A, f)
  k <- if (k_term) cc[length(cc)] else 0
  r <- if (k_term) cc[-length(cc)] else cc
  ## A zero residue would give an infinite precision for that block.
  r <- pmax(r, .Machine$double.eps * max(r))
  qq <- wl2_q(kind)
  b0 <- if (qq == 0) NULL else 1 / (1 + exp(-theta[1]))
  b <- exp(pmin(pmax(if (qq == 0) theta else theta[-1], -5), 60))
  o <- order(b)

  list(
    r = r[o], p = 1 - b[o], p0 = if (is.null(b0)) NULL else 1 - b0,
    q = qq, k = k,
    rel_err = sqrt(2 * best$cost / sum(cells$w * x^(2 * target))),
    x_min = x_min, theta = theta, kind = kind,
    alpha = alpha, d = d, m = m, type = type, s = s, k_term = k_term
  )
}

#' @name wl2_symbol
#' @title Evaluate the rational symbol of a weighted-L2 fit
#' @param cf Output of [rational.coefficients.wl2()].
#' @param x Points in (0, 1] at which to evaluate the symbol.
#' @return The values of the rational function.
#' @noRd
wl2_symbol <- function(cf, x) {
  s <- rowSums(outer(x, seq_along(cf$r), function(xx, i) {
    cf$r[i] * xx / (1 - cf$p[i] * xx)
  }))
  if (!is.null(cf$p0)) {
    q <- if (is.null(cf$q)) 1L else cf$q
    s <- (x / (1 - cf$p0 * x))^q * s
  }
  k <- if (is.null(cf$k)) 0 else cf$k
  s + k
}

#' @name wl2_alpha_grid
#' @title Grid of alpha values used for the weighted-L2 coefficient tables
#' @param m_alpha The integer part of alpha (covariance type only).
#' @param d Dimension.
#' @param type "covariance" or "operator".
#' @param by Resolution of the grid.
#' @return A vector of alpha values.
#' @noRd
wl2_alpha_grid <- function(m_alpha, d, type, by = 0.01) {
  alpha <- if (type == "operator") {
    seq(by, 2 - by, by = by)
  } else {
    m_alpha + seq(by, 1 - by, by = by)
  }
  ## alpha = nu + d / 2 with nu > 0, and the fit needs a finite trace
  alpha[alpha > d / 2 + by / 2]
}

#' @name wl2_check_table
#' @title Is a supplied table the right one for this model?
#' @description A table is fitted for one dimension, one order and one unit
#' interval of `alpha`, and using the wrong one would silently give a model
#' that is not the one asked for.
#' @param tab The table.
#' @param type,d,m,alpha What the model needs.
#' @return No return value; called for the error it raises.
#' @noRd
wl2_check_table <- function(tab, type, d, m, alpha) {
  if (!is.data.frame(tab) || is.null(attr(tab, "type"))) {
    stop("wl2_table is not a table produced by rspde.wl2.table().")
  }
  if (isTRUE(attr(tab, "k_term"))) {
    stop(paste0(
      "wl2_table was fitted with a constant term, which the covariance-based ",
      "construction has no block for. Such tables are for studying the ",
      "approximation, not for building a model."
    ))
  }
  got <- c(
    type = attr(tab, "type"), d = attr(tab, "d"), m = attr(tab, "m")
  )
  want <- c(type = type, d = d, m = m)
  bad <- names(want)[unlist(got[names(want)]) != unlist(want)]
  if (length(bad)) {
    stop(sprintf(
      "wl2_table was built for %s, but the model needs %s.",
      paste(sprintf("%s = %s", bad, unlist(got[bad])), collapse = ", "),
      paste(sprintf("%s = %s", bad, unlist(want[bad])), collapse = ", ")
    ))
  }
  rng <- range(tab$alpha)
  if (alpha < floor(rng[1]) || alpha > ceiling(rng[2])) {
    stop(sprintf(
      paste0(
        "wl2_table covers alpha in [%.2f, %.2f], but the model has ",
        "alpha = %.4f."
      ),
      rng[1], rng[2], alpha
    ))
  }
  invisible(NULL)
}

#' @name rspde.cache
#' @title Where generated coefficient tables are kept between sessions
#' @description The mesh-free weighted-L2 tables are shipped with the package,
#' but a table for a particular spectral interval (see [rspde.wl2.table()] and
#' the `kappa_ref` argument of [matern.operators()]) has to be computed, which
#' takes of the order of a second. Within a session such tables are cached in
#' memory. Enabling this cache stores them on disk as well, so that a second
#' session, or a script run again, finds them instead of refitting.
#'
#' The cache is off by default: a package should not write outside the session
#' temporary directory unless asked to. Calling `rspde.cache(TRUE)` is that
#' request, and creates the directory. The environment variable
#' `RSPDE_CACHE_DIR` sets it for non-interactive use, for instance from
#' `.Renviron`.
#'
#' Cached files are organised by package version, so upgrading the package does
#' not read tables produced by an older fit. Within a development version,
#' clear the cache after changing anything that affects the coefficients.
#'
#' @param dir `TRUE` to enable the cache at the standard location for this
#' package, `tools::R_user_dir("rSPDE", "cache")`; `FALSE` to disable it; a
#' path to use that directory instead; or missing to query the current setting.
#' @param clear Delete the tables that are currently stored?
#' @return The cache directory, or `NULL` if the cache is off, invisibly when
#' called for its effect.
#' @export
#' @seealso [rspde.wl2.table()], [matern.operators()]
#' @examples
#' rspde.cache()
#' \dontrun{
#' rspde.cache(TRUE) # store generated tables between sessions
#' rspde.cache(clear = TRUE) # throw away what is stored
#' rspde.cache(FALSE)
#' }
rspde.cache <- function(dir, clear = FALSE) {
  if (!missing(dir)) {
    if (isTRUE(dir)) {
      d <- tools::R_user_dir("rSPDE", "cache")
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      options(rSPDE.cache.dir = d)
    } else if (isFALSE(dir) || is.null(dir)) {
      options(rSPDE.cache.dir = NULL)
    } else {
      d <- path.expand(as.character(dir)[[1]])
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      options(rSPDE.cache.dir = d)
    }
  }
  d <- wl2_cache_dir()
  if (clear && !is.null(d)) {
    unlink(file.path(d, "wl2"), recursive = TRUE)
  }
  if (missing(dir) && !clear) {
    return(d)
  }
  invisible(d)
}

#' @name wl2_cache_dir
#' @title The directory tables are cached in, if any
#' @description The option set by [rspde.cache()] takes precedence, then the
#' environment variable `RSPDE_CACHE_DIR`. No default: the cache is off unless
#' the user turns it on.
#' @return A path, or `NULL`.
#' @noRd
wl2_cache_dir <- function() {
  d <- getOption("rSPDE.cache.dir", NULL)
  if (is.null(d)) {
    d <- Sys.getenv("RSPDE_CACHE_DIR", "")
    if (!nzchar(d)) {
      return(NULL)
    }
    d <- path.expand(d)
  }
  d
}

#' @name wl2_xmin_encode
#' @title Encode and decode `x_min` as a file name
#' @description A reversible encoding using only characters that are safe in a
#' path, so that the cached tables for one configuration can be listed and
#' their `x_min` read back without opening them.
#' @param x A value of `x_min`, or `NULL` for the mesh-free fit.
#' @return A string, or the value it encodes.
#' @noRd
wl2_xmin_encode <- function(x) {
  if (is.null(x)) {
    return("free")
  }
  e <- sprintf("%.12e", x)
  e <- gsub(".", "d", e, fixed = TRUE)
  e <- gsub("-", "m", e, fixed = TRUE)
  gsub("+", "p", e, fixed = TRUE)
}

#' @rdname wl2_xmin_encode
#' @noRd
wl2_xmin_decode <- function(s) {
  if (identical(s, "free")) {
    return(NULL)
  }
  e <- gsub("d", ".", s, fixed = TRUE)
  e <- gsub("m", "-", e, fixed = TRUE)
  e <- gsub("p", "+", e, fixed = TRUE)
  suppressWarnings(as.numeric(e))
}

#' @name wl2_cache_paths
#' @title Where a table is cached
#' @description One directory per configuration and one file per `x_min`
#' inside it, so that `wl2_cache_read()` can find a table fitted on a slightly
#' different mesh without opening every file. Cached files are separated by
#' package version, so that an upgrade never reads tables produced by an older
#' fit.
#' @param type,d,m,m_alpha,x_min,n_starts,n_grid,by,s,k_term The configuration.
#' @return A list with `dir` and `file`, or `NULL` if the cache is off.
#' @noRd
wl2_cache_paths <- function(type, d, m, m_alpha, x_min, n_starts, n_grid, by,
                            s = 0, k_term = FALSE) {
  root <- wl2_cache_dir()
  if (is.null(root)) {
    return(NULL)
  }
  dir <- file.path(
    root, "wl2", as.character(utils::packageVersion("rSPDE")),
    sprintf(
      "%s-d%s-m%s-a%s-n%s-g%s-b%s-w%s-k%s", type, d, m, m_alpha, n_starts,
      n_grid, gsub(".", "d", format(by), fixed = TRUE),
      gsub(".", "d", format(s), fixed = TRUE), if (k_term) 1 else 0
    )
  )
  list(dir = dir, file = file.path(dir, paste0(wl2_xmin_encode(x_min), ".rds")))
}

#' @name wl2_cache_match
#' @title A cached table that is valid for this spectral interval
#' @description The fit is optimal on \eqn{(x_{\min}, 1]}, so a table fitted
#' with a *smaller* `x_min` still covers the whole spectrum and remains a valid
#' model; one fitted with a larger `x_min` does not, and is never used. Among
#' the valid ones the largest `x_min` is the sharpest. Two meshes of the same
#' domain give values of `x_min` that differ by a fraction of a percent, and a
#' table half a decade too wide costs about a fifth of the error, so the search
#' stops at half the requested value; below that it is worth refitting.
#' @param dir The configuration's cache directory.
#' @param x_min The value the model needs.
#' @return A path, or `NULL`.
#' @noRd
wl2_cache_match <- function(dir, x_min) {
  if (is.null(x_min) || !dir.exists(dir)) {
    return(NULL)
  }
  files <- list.files(dir, pattern = "\\.rds$")
  if (length(files) == 0) {
    return(NULL)
  }
  xs <- vapply(sub("\\.rds$", "", files), function(f) {
    v <- wl2_xmin_decode(f)
    if (is.null(v) || length(v) != 1 || is.na(v)) NA_real_ else v
  }, numeric(1))
  ok <- which(!is.na(xs) & xs <= x_min & xs >= x_min / 2)
  if (length(ok) == 0) {
    return(NULL)
  }
  file.path(dir, files[ok[which.max(xs[ok])]])
}

#' @name wl2_cache_read
#' @title Read a cached table, if one is stored
#' @description Looks for a table fitted on exactly this interval, then for one
#' fitted on a slightly wider interval, which is still a valid model; see
#' `wl2_cache_match()`. Only the columns are stored, as for the tables shipped
#' with the package, so the attributes are reattached here, with the `x_min`
#' the table was actually fitted on rather than the one that was asked for.
#' @return The table, or `NULL`.
#' @noRd
wl2_cache_read <- function(type, d, m, m_alpha, x_min, n_starts, n_grid, by,
                           n_grid_quad, s = 0, k_term = FALSE) {
  paths <- wl2_cache_paths(type, d, m, m_alpha, x_min, n_starts, n_grid, by,
    s, k_term)
  if (is.null(paths)) {
    return(NULL)
  }
  f <- if (file.exists(paths$file)) {
    paths$file
  } else {
    wl2_cache_match(paths$dir, x_min)
  }
  if (is.null(f)) {
    return(NULL)
  }
  tab <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.data.frame(tab) || is.null(tab$rel_err)) {
    return(NULL)
  }
  got <- wl2_xmin_decode(sub("\\.rds$", "", basename(f)))
  wl2_table_attrs(tab, type, d, m, m_alpha, got, n_grid_quad, s, k_term)
}

#' @name wl2_cache_write
#' @title Store a table for later sessions
#' @description Failures are ignored: a cache that cannot be written is not a
#' reason for a model to fail.
#' @return No return value, called for the side effect.
#' @noRd
wl2_cache_write <- function(tab, type, d, m, m_alpha, x_min, n_starts, n_grid,
                            by, s = 0, k_term = FALSE) {
  paths <- wl2_cache_paths(type, d, m, m_alpha, x_min, n_starts, n_grid, by,
    s, k_term)
  if (is.null(paths)) {
    return(invisible(NULL))
  }
  bare <- tab
  attributes(bare) <- attributes(bare)[c("names", "row.names", "class")]
  try(
    {
      dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(bare, paths$file, compress = "xz")
    },
    silent = TRUE
  )
  invisible(NULL)
}

#' @name rspde.wl2.table
#' @title Weighted-L2 coefficients for a given spectral interval
#' @description Builds a table of weighted-L2 rational coefficients for a
#' specific lower end of the spectral interval, so that the approximation is
#' optimal on the interval the discrete operator actually has rather than on
#' \eqn{(0,1]}.
#'
#' The mesh-free tables that [matern.operators()] uses by default do not depend
#' on the mesh or on \eqn{\kappa} and are stored in the package, which is why
#' they cost nothing at set-up. They are, however, fitted on all of
#' \eqn{(0,1]}, while a discretised operator only ever sees
#' \eqn{x \in (x_{\min}, 1]} with \eqn{x_{\min}} set by the mesh resolution
#' and by \eqn{\kappa}. Fitting on that shorter interval spends the same
#' number of terms where they matter and is appreciably more accurate; the
#' price is that the table depends on the mesh and therefore has to be
#' computed, which takes of the order of a second.
#'
#' With the cache turned on (see [rspde.cache()]) there is usually nothing to
#' do by hand: a model asked for with `kappa_ref` computes its table once and
#' finds it again in later sessions without being told where it is. The two
#' routes share one store, so a table built here is what a later
#' `matern.operators(kappa_ref = ...)` picks up, and the other way round.
#'
#' This function is for what the cache does not cover: holding \eqn{x_{\min}}
#' fixed across several meshes or models rather than letting each derive its
#' own, inspecting the coefficients and the errors of the fit, or working
#' without writing to a cache at all. Pass the result to [matern.operators()]
#' as `wl2_table`; supplying it overrides `x_min` and `kappa_ref`. The result
#' is an ordinary data frame with attributes, so `saveRDS()` is also an option
#' if you would rather manage it yourself than enable the cache.
#'
#' Since \eqn{x_{\min}} grows with \eqn{\kappa}, a *lower bound* for
#' \eqn{\kappa} gives the shortest interval that is certainly long enough,
#' and hence a table that stays valid while \eqn{\kappa} is estimated. That
#' is what `kappa_ref` is.
#'
#' @param m The order of the rational approximation.
#' @param d The dimension of the domain. Taken from `mesh` when possible.
#' @param nu,alpha The smoothness. Give one of them; `alpha = nu + d/2`. Only
#' the integer part of `alpha` is used, since one table covers a whole unit
#' interval of `alpha`.
#' @param kappa_ref A lower bound for `kappa`, from which `x_min` is computed.
#' @param x_min The lower end of the spectral interval, if it is known; then
#' neither `kappa_ref` nor the mesh is needed.
#' @param C,G,mesh,loc_mesh The finite element matrices or the mesh, used with
#' `kappa_ref` to find `x_min`. See [rspde.xmin()].
#' @param type Either `"covariance"` or `"operator"`, matching the model the
#' table is for.
#' @param eigenvalue Passed to [rspde.xmin()].
#' @param s,k_term The weight exponent and the constant term of the fit; see
#' [rational.coefficients.wl2()]. These change what is being approximated and
#' are for studying the approximation rather than for building a model: a
#' constant term has no block in the covariance-based construction, and a
#' table carrying one is refused by [matern.operators()].
#' @param ... Passed to the fit, for instance `n_starts` or `by`.
#' @return A table of coefficients, with the attributes `type`, `d`, `m`,
#' `m_alpha`, `x_min` and `kind`, for use as the `wl2_table` argument of
#' [matern.operators()].
#' @export
#' @seealso [rspde.xmin()], [matern.operators()],
#' [rational.coefficients.wl2()]
#' @examples
#' \donttest{
#' mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 101))
#' tab <- rspde.wl2.table(m = 2, d = 1, nu = 0.8, kappa_ref = 5, mesh = mesh)
#' attr(tab, "x_min")
#' }
rspde.wl2.table <- function(m, d = NULL, nu = NULL, alpha = NULL,
                            kappa_ref = NULL, x_min = NULL,
                            C = NULL, G = NULL, mesh = NULL, loc_mesh = NULL,
                            type = c("covariance", "operator"),
                            eigenvalue = c("bound", "exact"),
                            s = 0, k_term = FALSE, ...) {
  type <- match.arg(type)
  if (is.null(d)) {
    d <- wl2_mesh_dim(mesh)
    if (is.null(d)) {
      stop("d could not be determined from the mesh; supply d.")
    }
  }
  if (is.null(alpha)) {
    if (is.null(nu)) {
      stop("Give either alpha or nu.")
    }
    alpha <- nu + d / 2
  }
  if (alpha <= 0) {
    stop("alpha must be positive.")
  }
  m_alpha <- floor(alpha)
  if (type == "covariance" && !(m_alpha %in% c(0, 1, 2))) {
    stop(paste0(
      "The weighted-L2 coefficients are only implemented for floor(alpha) ",
      "equal to 0, 1 or 2, but floor(alpha) = ", m_alpha, "."
    ))
  }
  if (is.null(x_min)) {
    if (is.null(kappa_ref) && is.null(mesh) && is.null(C)) {
      stop(paste0(
        "Give x_min, or kappa_ref together with a mesh (or C and G), so that ",
        "the spectral interval can be determined."
      ))
    }
    x_min <- rspde.xmin(
      C = C, G = G, mesh = mesh, kappa_ref = kappa_ref, nu = alpha - d / 2,
      loc_mesh = loc_mesh, eigenvalue = eigenvalue
    )
  }
  if (x_min <= 0 || x_min >= 1) {
    stop("x_min must lie in (0, 1).")
  }
  tab <- wl2_coefficient_table(
    d = d, m = m, m_alpha = if (type == "operator") floor(alpha) else m_alpha,
    x_min = x_min, type = type, s = s, k_term = k_term, ...
  )
  attr(tab, "kappa_ref") <- kappa_ref
  tab
}

#' @name wl2_mesh_dim
#' @title Dimension of the domain a mesh lives on
#' @param mesh An `fm_mesh_1d` or `fm_mesh_2d`, or `NULL`.
#' @return 1, 2, or `NULL` if it cannot be told.
#' @noRd
wl2_mesh_dim <- function(mesh) {
  if (is.null(mesh)) {
    return(NULL)
  }
  if (inherits(mesh, "fm_mesh_1d")) {
    return(1)
  }
  if (inherits(mesh, "fm_mesh_2d") || inherits(mesh, "inla.mesh")) {
    return(2)
  }
  NULL
}

#' @name wl2_table_attrs
#' @title Attach the descriptive attributes to a coefficient table
#' @description Everything about a table except its columns: what it is a table
#' of, and the quadrature the fits used. The quadrature is kept with the table
#' so that `wl2_interp_coefficients()` can re-solve the residues without
#' rebuilding it, and is reconstructed here rather than stored, since it
#' follows from `d`, `x_min` and `n_grid`.
#' @param tab The columns of the table.
#' @param type,d,m,m_alpha,x_min,n_grid,s,k_term The configuration it was built
#' for.
#' @return `tab`, with its attributes set.
#' @noRd
wl2_table_attrs <- function(tab, type, d, m, m_alpha, x_min, n_grid,
                            s = 0, k_term = FALSE) {
  attr(tab, "type") <- type
  attr(tab, "d") <- d
  attr(tab, "m") <- m
  attr(tab, "m_alpha") <- m_alpha
  attr(tab, "x_min") <- x_min
  attr(tab, "s") <- s
  attr(tab, "k_term") <- k_term
  attr(tab, "kind") <- if (type == "covariance") {
    c("plain", "shifted", "shifted2")[m_alpha + 1]
  } else {
    "plain"
  }
  cells <- wl2_weyl_cells(d, if (is.null(x_min)) 1e-26 else x_min,
    if (is.null(x_min)) max(n_grid, 900) else n_grid, s = s
  )
  attr(tab, "quad") <- list(x = cells$x, w = cells$w, sw = sqrt(cells$w))
  tab
}

#' @name wl2_shipped_table
#' @title A precomputed mesh-free coefficient table, if there is one
#' @description The mesh-free covariance tables do not depend on the mesh or on
#' `kappa`, and are the default for `type = "covariance"`, so they are built
#' once by `data-raw/wl2_tables.R` and shipped with the package. Rebuilding one
#' takes up to four seconds; looking it up costs the few milliseconds of the
#' quadrature. Tables for other configurations, or for non-default settings of
#' the fit, are not stored and are computed on demand.
#' @param d,m,m_alpha,n_grid The configuration.
#' @return The table, or `NULL` if none is stored for this configuration.
#' @noRd
wl2_shipped_table <- function(d, m, m_alpha, n_grid) {
  ## get0(), not the object itself: data-raw/wl2_tables.R runs against a
  ## namespace in which it does not exist yet.
  tabs <- get0("wl2_meshfree_tables",
    envir = environment(wl2_shipped_table), ifnotfound = NULL
  )
  tab <- tabs[[paste("covariance", d, m, m_alpha, sep = "_")]]
  if (is.null(tab)) {
    return(NULL)
  }
  wl2_table_attrs(tab, "covariance", d, m, m_alpha, NULL, n_grid, 0, FALSE)
}

#' @name wl2_theta_canonical
#' @title Outer parameters of a fit, in sorted order
#' @description The poles enter the basis symmetrically, so a fit's parameters
#' are only defined up to a permutation. Sorting them gives a representation
#' that varies smoothly along the continuation in `alpha`, which is what makes
#' it possible to extrapolate from one grid point to the next.
#' @param cf A fit, as returned by [rational.coefficients.wl2()].
#' @return The outer parameters: `logit(b_0)` (shifted classes only) followed by
#' `log(b)` in increasing order.
#' @noRd
wl2_theta_canonical <- function(cf) {
  v <- log(1 - cf$p)
  if (!is.null(cf$p0)) {
    v <- c(stats::qlogis(1 - cf$p0), v)
  }
  pmin(pmax(v, -5), 60)
}

#' @name wl2_restart_starts
#' @title Warm restarts for a continuation step that went wrong
#' @description A fit whose smallest residue has been driven to zero is really
#' a solution of order `m - 1`. The natural repair is the move the continuation
#' in `m` already uses: drop the dead term and put a pole back into each gap of
#' the remaining ones. The poles enter the basis symmetrically, so they can be
#' taken in the sorted order reported by the fit.
#' @param cf A fit, as returned by [rational.coefficients.wl2()].
#' @param shifted Does the class have a shift, that is, is `q` positive?
#' @return A list of starting values for the outer parameters.
#' @noRd
wl2_restart_starts <- function(cf, shifted) {
  lb <- log(1 - cf$p)[-which.min(cf$r)]
  lb <- pmin(pmax(lb, -5), 60)
  wl2_insert_pole(
    if (shifted) c(stats::qlogis(1 - cf$p0), lb) else lb, shifted
  )
}

#' @name wl2_coefficient_table
#' @title Table of weighted-L2 rational coefficients
#' @description Computes the coefficients on a grid of `alpha` values, using
#' continuation in `alpha` (the fit for one grid point is used as starting
#' value for the next) so that a full multistart is only needed for the first
#' point. Results are cached, so that repeatedly creating models with the same
#' configuration is cheap.
#' @param d Dimension of the domain.
#' @param m Order of the rational approximation.
#' @param m_alpha The integer part of alpha. Ignored for `type = "operator"`.
#' @param x_min Lower end of the spectral interval, or `NULL` for the mesh-free
#' fit.
#' @param type "covariance" or "operator".
#' @param n_starts Number of starts used for the first grid point.
#' @param n_grid Number of quadrature nodes.
#' @param by Resolution of the alpha grid.
#' @param cache Should the result be cached?
#' @param s,k_term The weight exponent and the constant term of the fit; see
#' [rational.coefficients.wl2()]. Tables with either set are never taken from
#' the package's stored tables, which are for `s = 0` and no constant term.
#' @param shipped May an already computed table be used: one stored in the
#' package, where there is one for exactly this configuration, or one written
#' to the persistent cache by an earlier session (see [rspde.cache()])?
#' `FALSE` forces the fits to be run, which is what `data-raw/wl2_tables.R` and
#' the test that checks the stored tables do.
#' @return A data frame with columns `alpha`, `r1`, ..., `p1`, ..., `k`, and
#' `p0` for the shifted classes, together with the attributes `type`, `d`, `m`,
#' `m_alpha`, `x_min` and `kind`.
#' @noRd
wl2_coefficient_table <- function(d, m, m_alpha, x_min = NULL,
                                  type = "covariance", n_starts = 8,
                                  n_grid = 300, by = 0.01, cache = TRUE,
                                  shipped = TRUE, s = 0, k_term = FALSE) {
  k_term <- isTRUE(k_term)
  key <- paste(type, d, m, m_alpha, if (is.null(x_min)) "free" else
    signif(x_min, 12), n_starts, n_grid, by, s, k_term, sep = "_")
  if (cache && !is.null(wl2_cache[[key]])) {
    return(wl2_cache[[key]])
  }

  ## A stored table, where there is one for exactly this fit.
  ## The stored tables are for the default weight and no constant term.
  if (shipped && is.null(x_min) && type == "covariance" &&
    n_starts == 8 && n_grid == 300 && by == 0.01 && s == 0 && !k_term) {
    tab <- wl2_shipped_table(d, m, m_alpha, n_grid)
    if (!is.null(tab)) {
      if (cache) {
        wl2_cache[[key]] <- tab
      }
      return(tab)
    }
  }

  ## A table this or an earlier session generated, if the user has turned the
  ## persistent cache on. This is what makes a mesh-specific table a one-off
  ## rather than a per-session cost.
  if (cache && shipped) {
    tab <- wl2_cache_read(
      type, d, m, m_alpha, x_min, n_starts, n_grid, by, n_grid, s, k_term
    )
    if (!is.null(tab)) {
      wl2_cache[[key]] <- tab
      return(tab)
    }
  }

  alphas <- wl2_alpha_grid(m_alpha, d, type, by)
  if (length(alphas) == 0) {
    stop("No valid alpha values for the requested weighted-L2 table.")
  }
  n_col <- if (type == "operator") m + 1 else m

  r_mat <- matrix(NA_real_, length(alphas), n_col)
  p_mat <- matrix(NA_real_, length(alphas), n_col)
  p0_vec <- rep(NA_real_, length(alphas))
  k_vec <- rep(0, length(alphas))
  err_vec <- rep(NA_real_, length(alphas))
  theta_list <- vector("list", length(alphas))

  kind_shifted <- type == "covariance" && m_alpha > 0
  is_collapsed <- function(cf) min(cf$r) < 1e-8 * max(cf$r)

  prev <- NULL
  ## Was the collapse at the previous grid point confirmed to be the optimum?
  ## Near an end of the alpha range the best fit of order m genuinely has fewer
  ## than m active terms, and then every repair returns the same solution;
  ## re-deriving that at each of the following grid points is the single most
  ## expensive thing the sweep can do.
  settled <- FALSE
  ## The optima move smoothly with alpha, so the previous two give a linearly
  ## extrapolated start that is centred on the new point rather than one grid
  ## step behind it. Interpolating a coarse first pass onto the intermediate
  ## points is equivalent and costs the same; both save about a tenth of the
  ## sweep, which is all there is to save once the start is this good.
  prev2 <- NULL
  for (i in seq_along(alphas)) {
    start_i <- if (is.null(prev)) {
      NULL
    } else if (is.null(prev2)) {
      wl2_theta_canonical(prev)
    } else {
      2 * wl2_theta_canonical(prev) - wl2_theta_canonical(prev2)
    }
    cf <- rational.coefficients.wl2(
      alpha = alphas[i], d = d, m = m, x_min = x_min, type = type,
      n_starts = if (is.null(prev)) n_starts else 1L,
      n_grid = n_grid, start = start_i, s = s, k_term = k_term
    )
    ## Guard against the continuation ending up in a worse local minimum, or
    ## collapsing to the optimum of a smaller order (a residue driven to zero,
    ## which wastes a block and then propagates through the continuation).
    collapsed <- is_collapsed(cf)
    worse <- !is.null(prev) && cf$rel_err > 2 * prev$rel_err
    if (!is.null(prev) && (worse || (collapsed && !settled))) {
      ## Re-inserting the dead pole is some five times cheaper than a cold
      ## multistart and has found the same optimum in every case tried.
      best <- NULL
      for (st in wl2_restart_starts(cf, kind_shifted)) {
        cf2 <- tryCatch(
          rational.coefficients.wl2(
            alpha = alphas[i], d = d, m = m, x_min = x_min, type = type,
            n_starts = 1L, n_grid = n_grid, start = st, s = s, k_term = k_term
          ),
          error = function(e) NULL
        )
        if (!is.null(cf2) && (is.null(best) || cf2$rel_err < best$rel_err)) {
          best <- cf2
        }
      }
      ## A cold multistart is kept as the last resort, but only for a point
      ## that is actually anomalous: a collapse that the restarts cannot undo
      ## is a property of the problem, not a failure of the continuation.
      if (worse && (is.null(best) || best$rel_err >= cf$rel_err)) {
        best <- rational.coefficients.wl2(
          alpha = alphas[i], d = d, m = m, x_min = x_min, type = type,
          n_starts = n_starts, n_grid = n_grid, s = s, k_term = k_term
        )
      }
      if (!is.null(best) && (best$rel_err < cf$rel_err ||
        (collapsed && !is_collapsed(best)))) {
        cf <- best
      }
      settled <- is_collapsed(cf)
    } else {
      settled <- collapsed && settled
    }
    r_mat[i, ] <- cf$r
    p_mat[i, ] <- cf$p
    k_vec[i] <- cf$k
    if (!is.null(cf$p0)) {
      p0_vec[i] <- cf$p0
    }
    err_vec[i] <- cf$rel_err
    theta_list[[i]] <- cf$theta
    prev2 <- prev
    prev <- cf
  }

  ## Polish: with a monotone trend in alpha the error at a grid point cannot
  ## exceed both of its neighbours, so a point that does is a local minimum of
  ## the outer problem. Refit it warm-started from either neighbour.
  for (pass in 1:2) {
    interior <- seq_along(alphas)[-c(1, length(alphas))]
    spikes <- interior[err_vec[interior] >
      1.5 * pmax(err_vec[interior - 1], err_vec[interior + 1])]
    if (length(spikes) == 0) {
      break
    }
    for (i in spikes) {
      for (j in c(i - 1, i + 1)) {
        cf <- rational.coefficients.wl2(
          alpha = alphas[i], d = d, m = m, x_min = x_min, type = type,
          n_starts = 1L, n_grid = n_grid, start = theta_list[[j]],
          s = s, k_term = k_term
        )
        if (cf$rel_err < err_vec[i]) {
          r_mat[i, ] <- cf$r
          p_mat[i, ] <- cf$p
          k_vec[i] <- cf$k
          if (!is.null(cf$p0)) {
            p0_vec[i] <- cf$p0
          }
          err_vec[i] <- cf$rel_err
          theta_list[[i]] <- cf$theta
        }
      }
    }
  }

  tab <- data.frame(alpha = alphas)
  for (j in seq_len(n_col)) {
    tab[[paste0("r", j)]] <- r_mat[, j]
  }
  for (j in seq_len(n_col)) {
    tab[[paste0("p", j)]] <- p_mat[, j]
  }
  tab$k <- k_vec
  if (all(!is.na(p0_vec))) {
    tab$p0 <- p0_vec
  }
  tab$rel_err <- err_vec
  tab <- wl2_table_attrs(tab, type, d, m, m_alpha, x_min, n_grid, s, k_term)
  if (cache) {
    wl2_cache[[key]] <- tab
    wl2_cache_write(tab, type, d, m, m_alpha, x_min, n_starts, n_grid, by,
      s, k_term)
  }
  tab
}

## Cache of the coefficient tables, keyed by the configuration.
wl2_cache <- new.env(parent = emptyenv())

#' @name wl2_clear_cache
#' @title Empty the cache of weighted-L2 coefficient tables
#' @return No return value, called for the side effect.
#' @noRd
wl2_clear_cache <- function() {
  rm(list = ls(wl2_cache), envir = wl2_cache)
  invisible(NULL)
}

#' @name wl2_interp_coefficients
#' @title Interpolate the weighted-L2 coefficients
#' @description Interpolates a table produced by `wl2_coefficient_table()` at a
#' given value of `alpha`. The poles are obtained by spline interpolation of the
#' internal parameters (`log(b)`, and `logit(b_0)` for the shifted class), which
#' are the smooth ones, and the residues are then re-solved exactly by
#' non-negative least squares for the interpolated poles. This is as accurate as
#' a direct fit, also close to integer `alpha` where the residues and poles of
#' the table vary over many orders of magnitude.
#' @param tab The table.
#' @param alpha The value of alpha.
#' @return A list with elements `r`, `p`, `p0`, `k` and `rel_err`.
#' @noRd
wl2_interp_coefficients <- function(tab, alpha) {
  alpha <- min(max(alpha, min(tab$alpha)), max(tab$alpha))
  n_col <- sum(grepl("^r[0-9]+$", names(tab)))
  kind <- attr(tab, "kind")
  sp <- function(y) stats::spline(tab$alpha, y, xout = alpha)$y

  theta <- vapply(seq_len(n_col), function(i) {
    sp(log(1 - tab[[paste0("p", i)]]))
  }, numeric(1))
  qq <- wl2_q(kind)
  if (qq > 0) {
    b0 <- pmin(pmax(1 - tab$p0, 1e-12), 1 - 1e-12)
    theta <- c(sp(log(b0 / (1 - b0))), theta)
  }

  quad <- attr(tab, "quad")
  k_term <- isTRUE(attr(tab, "k_term"))
  target <- if (identical(attr(tab, "type"), "operator")) alpha / 2 else alpha
  f <- quad$x^target * quad$sw
  A <- wl2_basis(kind, quad$x, theta, k_term) * quad$sw
  ## The residues, and the constant term with them, are re-solved rather than
  ## interpolated: given the poles this is the exact optimum, and it keeps the
  ## interpolated coefficients a genuine non-negative fit.
  cc <- wl2_nnls(A, f)
  res <- as.vector(A %*% cc) - f
  k <- if (k_term) cc[length(cc)] else 0
  r <- if (k_term) cc[-length(cc)] else cc
  ## A zero residue would give an infinite precision for that block.
  r <- pmax(r, .Machine$double.eps * max(r))

  b <- exp(pmin(pmax(if (qq > 0) theta[-1] else theta, -5), 60))
  list(
    r = r, p = 1 - b,
    p0 = if (qq > 0) 1 - 1 / (1 + exp(-theta[1])) else NULL,
    q = qq, k = k, theta = theta,
    rel_err = sqrt(sum(res^2) / sum(quad$w * quad$x^(2 * target)))
  )
}

#' @name wl2_domain_diameter
#' @title Diameter of the domain of a mesh
#' @description The diagonal of the bounding box of the mesh nodes.
#' @param mesh An `fmesher` mesh, a `metric_graph` object, or `NULL`.
#' @param loc_mesh Mesh locations, used if `mesh` is `NULL`.
#' @return The diameter, or `NULL` if it cannot be determined.
#' @noRd
wl2_domain_diameter <- function(mesh = NULL, loc_mesh = NULL) {
  loc <- NULL
  if (!is.null(mesh)) {
    if (inherits(mesh, "metric_graph")) {
      loc <- mesh$get_vertices()
    } else if (inherits(mesh, "fm_mesh_1d")) {
      loc <- matrix(mesh$loc, ncol = 1)
    } else if (!is.null(mesh$loc)) {
      loc <- as.matrix(mesh$loc)
    }
  } else if (!is.null(loc_mesh)) {
    loc <- as.matrix(loc_mesh)
  }
  if (is.null(loc) || length(loc) == 0) {
    return(NULL)
  }
  sqrt(sum((apply(loc, 2, max) - apply(loc, 2, min))^2))
}

#' @name wl2_lambda_max
#' @title Largest eigenvalue of a sparse symmetric matrix
#' @description Lanczos with full reorthogonalisation. How fast this converges
#' depends on how separated the top of the spectrum is. Ritz values are interior 
#' to the spectrum, so this converges from *below*, which is the unsafe direction: 
#' too small a \eqn{\lambda_{\max}} gives too large an \eqn{x_{\min}} and leaves 
#' the top of the spectrum outside the fitted interval. The caller therefore 
#' allows a one per cent margin and caps the result at the Gershgorin bound, so 
#' that the value used is never below the true largest eigenvalue and never more 
#' than one per cent above it.
#' @param S A sparse symmetric matrix.
#' @param k Number of Lanczos steps.
#' @return The largest eigenvalue, or `NA_real_` if the iteration breaks down.
#' @noRd
wl2_lambda_max <- function(S, k = 20) {
  n <- dim(S)[1]
  if (n <= k + 2) {
    return(max(eigen(as.matrix((S + Matrix::t(S)) / 2),
      symmetric = TRUE, only.values = TRUE
    )$values))
  }
  ## A fixed starting vector, so that x_min does not depend on the state of the
  ## user's random number generator.
  q <- rep_len(c(1, -1, 1, 1, -1), n) + seq_len(n) / n
  q <- q / sqrt(sum(q^2))
  Q <- matrix(0, n, k)
  al <- be <- numeric(k)
  qm <- numeric(n)
  b <- 0
  kk <- k
  for (j in seq_len(k)) {
    Q[, j] <- q
    v <- as.vector(S %*% q)
    al[j] <- sum(q * v)
    v <- v - al[j] * q - b * qm
    ## Full reorthogonalisation: k is small, and without it the Ritz values
    ## acquire spurious copies of the dominant eigenvalue. Projecting with the
    ## whole of Q rather than its first j columns gives the same result, since
    ## the remaining columns are still zero. It does more arithmetic but copies
    ## nothing, and copying is what costs: 2.3 times faster at 10^4 nodes and
    ## above, against about 25% slower below 5000, which is the right trade for
    ## a routine whose absolute cost only matters on large meshes.
    v <- v - Q %*% crossprod(Q, v)
    b <- sqrt(sum(v^2))
    if (!is.finite(b) || b < 1e-12) {
      kk <- j
      break
    }
    be[j] <- b
    qm <- q
    q <- as.vector(v) / b
  }
  Tm <- diag(al[seq_len(kk)], kk)
  if (kk > 1) {
    i <- seq_len(kk - 1)
    Tm[cbind(i, i + 1)] <- be[i]
    Tm[cbind(i + 1, i)] <- be[i]
  }
  ev <- tryCatch(max(eigen(Tm, symmetric = TRUE, only.values = TRUE)$values),
    error = function(e) NA_real_
  )
  if (length(ev) != 1 || !is.finite(ev)) NA_real_ else ev
}


#' Lower end of the spectral interval of a rational approximation
#'
#' Computes the value \eqn{x_{\min} = 1/(1 + \mu_{\max}/\kappa_{\mathrm{lo}}^2)}
#' used by the weighted-L2 rational coefficients, where \eqn{\mu_{\max}} is the
#' largest generalised eigenvalue of \eqn{(G, C)} and
#' \eqn{\kappa_{\mathrm{lo}}} is a lower bound for the values that \eqn{\kappa}
#' can take.
#'
#' The reference \eqn{\kappa_{\mathrm{lo}}} must be a *lower* bound: a reference
#' above the true \eqn{\kappa} leaves part of the spectrum outside the fitted
#' interval, and the error then grows by one to two orders of magnitude. A
#' reference below the true \eqn{\kappa} is safe; the approximation then
#' degrades gracefully towards the mesh-free fit. The default,
#' \eqn{\kappa_{\mathrm{lo}} = \sqrt{8\nu}/\mathrm{diam}}, corresponds to a
#' range equal to the diameter of the domain.
#'
#' @param C The mass matrix of the finite element discretisation.
#' @param G The stiffness matrix of the finite element discretisation.
#' @param mesh An optional mesh; `C` and `G` are computed from it if they are
#' not given, and it is used for the diameter of the domain.
#' @param kappa_ref The reference value \eqn{\kappa_{\mathrm{lo}}}. If `NULL`,
#' it is taken to be `sqrt(8 * nu) / diameter`.
#' @param nu The smoothness parameter, used for the default `kappa_ref`.
#' @param diameter The diameter of the domain, used for the default
#' `kappa_ref`. Computed from `mesh` if not given.
#' @param loc_mesh Mesh locations, an alternative to `mesh` for the diameter.
#' @param eigenvalue Either `"bound"` (the default), which uses the Gershgorin
#' bound `max(rowSums(abs(G)) / diag(C))` for \eqn{\mu_{\max}}, or `"exact"`,
#' which computes the eigenvalue by Lanczos iteration. The bound is an
#' over-estimate of \eqn{\mu_{\max}}, and therefore gives a conservative (too
#' small) \eqn{x_{\min}}, which is the safe direction.
#' @return The value of \eqn{x_{\min}}.
#' @export
#' @seealso [rational.coefficients.wl2()]
#' @examples
#' x <- seq(from = 0, to = 1, length.out = 201)
#' fem <- rSPDE.fem1d(x)
#' # reference kappa corresponding to a range equal to the domain diameter
#' rspde.xmin(C = fem$C, G = fem$G, loc_mesh = x, nu = 0.5)
#' # a user-supplied lower bound for kappa
#' rspde.xmin(C = fem$C, G = fem$G, kappa_ref = 10)
rspde.xmin <- function(C = NULL, G = NULL, mesh = NULL, kappa_ref = NULL,
                       nu = NULL, diameter = NULL, loc_mesh = NULL,
                       eigenvalue = c("bound", "exact")) {
  eigenvalue <- match.arg(eigenvalue)
  if ((is.null(C) || is.null(G))) {
    if (is.null(mesh)) {
      stop("Either C and G, or mesh, must be provided.")
    }
    fem <- fm_fem(mesh)
    C <- fem$c0
    G <- fem$g1
  }
  cdiag <- Matrix::rowSums(C)
  if (is.null(kappa_ref)) {
    if (is.null(diameter)) {
      diameter <- wl2_domain_diameter(mesh, loc_mesh)
    }
    if (is.null(diameter) || is.null(nu) || diameter <= 0) {
      stop(paste0(
        "kappa_ref could not be determined; supply kappa_ref, or a mesh ",
        "(or loc_mesh or diameter) together with nu."
      ))
    }
    kappa_ref <- sqrt(8 * max(nu, 1e-6)) / diameter
  }
  kappa_ref <- min(kappa_ref)
  if (kappa_ref <= 0) {
    stop("kappa_ref must be positive.")
  }
  ## Gershgorin: always an upper bound for the largest eigenvalue, so it always
  ## gives an x_min that is safe, if in two dimensions about 70% too small.
  mu_bound <- max(Matrix::rowSums(abs(G)) / cdiag)
  mu_max <- mu_bound
  if (eigenvalue == "exact") {
    Ci2 <- Matrix::Diagonal(length(cdiag), 1 / sqrt(cdiag))
    mu_lanczos <- wl2_lambda_max(Ci2 %*% G %*% Ci2)
    if (is.finite(mu_lanczos) && mu_lanczos > 0) {
      ## Lanczos converges from below; the margin and the cap keep the result
      ## between the true eigenvalue and one per cent above it.
      mu_max <- min(mu_bound, 1.01 * mu_lanczos)
    } else {
      warning(paste0(
        "The largest eigenvalue could not be computed; using the Gershgorin ",
        "bound instead."
      ))
    }
  }
  1 / (1 + mu_max / kappa_ref^2)
}

#' @name wl2_polymul
#' @title Multiply two polynomials given by their coefficients
#' @param a,b Coefficient vectors in ascending order.
#' @return The coefficients of the product, in ascending order.
#' @noRd
wl2_polymul <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1)
  for (i in seq_along(a)) {
    out[i + seq_along(b) - 1] <- out[i + seq_along(b) - 1] + a[i] * b
  }
  out
}

#' @name wl2_roots
#' @title Convert weighted-L2 partial fractions to the operator-based form
#' @description Converts \eqn{\sum_i s_i/(\lambda - \lambda_i)} into the form
#' `factor` \eqn{\prod_i (1 - rc_i \lambda) / \prod_j (1 - rb_j \lambda)} used by
#' [fractional.operators()].
#' @param cf Output of [rational.coefficients.wl2()] or
#' `wl2_interp_coefficients()`, with `m + 1` terms.
#' @param tol Tolerance used when checking that the two forms agree.
#' @return A list with elements `rb`, `rc` and `factor`, as returned by
#' `get.roots()`.
#' @noRd
wl2_roots <- function(cf, tol = 1e-8) {
  if (!is.null(cf$k) && cf$k != 0) {
    stop(paste0(
      "The operator factorisation assumes no constant term, but these ",
      "coefficients have k = ", signif(cf$k, 4), ". Coefficients fitted with ",
      "k_term = TRUE are for evaluating the symbol, not for building an ",
      "operator-based model."
    ))
  }
  lambda <- cf$p
  s <- cf$r
  M <- length(lambda)
  ## N(lambda) = sum_i s_i prod_{j != i} (lambda - lambda_j), of degree M - 1
  Ncoef <- numeric(M)
  for (i in seq_len(M)) {
    pc <- 1
    for (j in seq_len(M)[-i]) {
      pc <- wl2_polymul(pc, c(-lambda[j], 1))
    }
    Ncoef <- Ncoef + s[i] * pc
  }
  rts <- polyroot(Ncoef)
  if (max(abs(Im(rts))) > 1e-6 * max(1, max(abs(Re(rts))))) {
    stop("The numerator of the rational approximation has complex roots.")
  }
  rts <- sort(Re(rts), decreasing = TRUE)
  if (any(rts >= 0)) {
    stop("The numerator of the rational approximation has non-negative roots.")
  }

  rb <- 1 / lambda
  rc <- 1 / rts
  factor <- Ncoef[1] / prod(-lambda)

  ## The two forms must agree; check on a grid of the spectral interval.
  lam <- exp(seq(0, log(1 / max(cf$x_min, 1e-20)), length.out = 50))
  lhs <- factor * apply(outer(lam, rc, function(l, r) 1 - r * l), 1, prod) /
    apply(outer(lam, rb, function(l, r) 1 - r * l), 1, prod)
  rhs <- rowSums(outer(lam, seq_along(s), function(l, i) s[i] / (l - lambda[i])))
  if (max(abs(lhs - rhs)) > tol * max(abs(rhs))) {
    stop("The conversion to the operator-based form failed.")
  }

  list(rb = rb, rc = rc, factor = factor)
}

#' @name wl2_diag_inverse
#' @title Diagonal of the inverse of a sparse SPD matrix
#' @description Uses `INLA::inla.qinv()` when INLA is available, and otherwise
#' a direct solve, which is only feasible for moderate sizes.
#' @param Q A sparse symmetric positive definite matrix.
#' @param max_direct Largest dimension for which the direct solve is attempted.
#' @return The diagonal of the inverse.
#' @noRd
wl2_diag_inverse <- function(Q, max_direct = 4000) {
  n <- dim(Q)[1]
  if (requireNamespace("INLA", quietly = TRUE)) {
    return(Matrix::diag(INLA::inla.qinv(Matrix::forceSymmetric(Q))))
  }
  if (n > max_direct) {
    stop(paste0(
      "The nodal variance correction needs the diagonal of the inverse of a ",
      n, " x ", n, " matrix. Install INLA for a sparse (Takahashi) ",
      "computation, or use variance_correction = 'none'."
    ))
  }
  Matrix::diag(Matrix::solve(Matrix::forceSymmetric(Q), Matrix::Diagonal(n)))
}

#' @name wl2_nodal_correction
#' @title Nodal variance correction for weighted-L2 covariance models
#' @description Without a constant term the nodal variance of the
#' approximation is not exactly \eqn{\sigma^2}; this returns the diagonal
#' \eqn{D = \sigma^2 - \mathrm{diag}(\Sigma_{\mathrm{approx}})} that restores
#' it. In observation models, adding `D` to the covariance is equivalent to
#' adding `A D A^T` to the covariance of the measurement noise.
#' @param Q The precision matrix of the latent blocks.
#' @param n_blocks The number of blocks.
#' @param sigma The target marginal standard deviation.
#' @return The correction vector, of length `nrow(Q) / n_blocks`.
#' @noRd
wl2_nodal_correction <- function(Q, n_blocks, sigma) {
  if (is.list(Q)) {
    dg <- Reduce(`+`, lapply(Q, wl2_diag_inverse))
    return(sigma^2 - dg)
  }
  n <- dim(Q)[1] / n_blocks
  dg <- numeric(n)
  for (j in seq_len(n_blocks)) {
    idx <- (j - 1) * n + seq_len(n)
    dg <- dg + wl2_diag_inverse(Q[idx, idx, drop = FALSE])
  }
  sigma^2 - dg
}

#' Recompute the weighted-L2 coefficients for an estimated kappa
#'
#' Two-stage use of the weighted-L2 rational coefficients. The coefficients
#' must not depend on `kappa` during estimation, since that would make the
#' likelihood non-smooth in `kappa`, so they are computed at set-up for a
#' conservative reference `kappa`. Once `kappa` has been estimated, this
#' function recomputes them once for a reference derived from the estimate, for
#' the final likelihood evaluation and for prediction.
#'
#' The reference must remain a lower bound for `kappa`, which is why the
#' default divides the estimate by `safety`. A reference above the true `kappa`
#' leaves part of the spectrum outside the fitted interval and the error grows
#' by one to two orders of magnitude, whereas a reference below it degrades
#' gracefully towards the mesh-free fit.
#'
#' @param object A model created by [matern.operators()] with
#' `type_rational_approximation = "wl2"`.
#' @param kappa_ref The new reference value. If `NULL`, the `kappa` of the
#' object divided by `safety` is used.
#' @param safety The safety factor applied to the `kappa` of the object when
#' `kappa_ref` is not given.
#' @param ... Further arguments passed to `update()`.
#' @return The updated model.
#' @export
#' @seealso [rational.coefficients.wl2()], [rspde.xmin()]
#' @examples
#' x <- seq(from = 0, to = 1, length.out = 101)
#' op <- matern.operators(
#'   loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
#'   parameterization = "matern", type = "operator",
#'   type_rational_approximation = "wl2"
#' )
#' op <- update_rational_coefficients(op)
update_rational_coefficients <- function(object, kappa_ref = NULL,
                                         safety = 3, ...) {
  if (!identical(object$type_rational_approximation, "wl2")) {
    stop(paste0(
      "update_rational_coefficients() is only meaningful for models with ",
      "type_rational_approximation = 'wl2'."
    ))
  }
  if (is.null(kappa_ref)) {
    if (is.null(object$kappa)) {
      stop("The object does not contain kappa; supply kappa_ref.")
    }
    kappa_ref <- min(object$kappa) / safety
  }
  stats::update(object, kappa_ref = kappa_ref, x_min = NULL, ...)
}

#' @name wl2_cgeneric_table
#' @title A dense weighted-L2 table for the INLA cgeneric interface
#' @description The cgeneric model looks its coefficients up in a table rather
#' than fitting anything, so the weighted-L2 coefficients reach it the same way
#' the tabulated ones do: evaluated in R on a fine grid of `alpha` and passed
#' as a matrix.
#'
#' The tabulated types need only 999 rows, indexed by the fractional part of
#' `alpha`, because they approximate \eqn{x^{\{\alpha\}}} and the integer
#' factor \eqn{x^{\lfloor\alpha\rfloor}} is exact and built separately. The
#' weighted-L2 classes are fitted to \eqn{x^\alpha} as a whole, with the
#' integer factor carrying the shift \eqn{p_0}, so their coefficients depend on
#' all of `alpha`. The table therefore holds one block of 999 rows per
#' \eqn{\lfloor\alpha\rfloor} reachable from the prior on `nu`, and the row for
#' a given `alpha` is
#' `(floor(alpha) - m_alpha_min) * 999 + round(1000 * frac(alpha))`.
#'
#' The columns are those of `get_rational_coefficients()` with the constant
#' term replaced by the shift: `alpha`, `r1..rm`, `p1..pm`, `p0`. Each row
#' costs one interpolation, about 0.3 ms, so a block takes a third of a second.
#' @param d The dimension of the domain.
#' @param m The order of the rational approximation.
#' @param nu_upper_bound The upper bound of the prior on `nu`.
#' @param wl2_table An optional mesh-free table per band; by default the
#' shipped one.
#' @return A matrix, with the attribute `m_alpha_min`.
#' @noRd
wl2_cgeneric_table <- function(d, m, nu_upper_bound, wl2_table = NULL) {
  ## alpha = nu + d / 2 with nu in (0, nu_upper_bound)
  m_alpha_min <- floor(d / 2 + 1e-10)
  ## nu is strictly below its upper bound, so alpha is strictly below
  ## nu_upper_bound + d/2 and the top band is the one just under it.
  m_alpha_max <- floor(nu_upper_bound + d / 2 - 1e-10)
  if (m_alpha_max > 2) {
    stop(paste0(
      "type.rational.approx = 'wl2' is only implemented for floor(alpha) ",
      "equal to 0, 1 or 2, but nu.upper.bound = ", nu_upper_bound,
      " with d = ", d, " reaches floor(alpha) = ", m_alpha_max, "."
    ))
  }
  frac <- seq_len(999) / 1000
  out <- NULL
  for (ma in m_alpha_min:m_alpha_max) {
    tab <- if (is.null(wl2_table)) {
      wl2_coefficient_table(d = d, m = m, m_alpha = ma, type = "covariance")
    } else {
      wl2_table
    }
    rows <- vapply(ma + frac, function(a) {
      cf <- wl2_interp_coefficients(tab, a)
      c(a, cf$r, cf$p, if (is.null(cf$p0)) 0 else cf$p0)
    }, numeric(2 * m + 2))
    out <- rbind(out, t(rows))
  }
  colnames(out) <- c(
    "alpha", paste0("r", seq_len(m)), paste0("p", seq_len(m)), "p0"
  )
  attr(out, "m_alpha_min") <- m_alpha_min
  out
}

#' @name wl2_operator_type
#' @title The rational type to use for the operator-based construction
#' @description The operator-based factorisation has a single table of roots,
#' produced by the chebfun lower-bound method, so `"brasil"` and `"chebfun"`
#' are refused there rather than quietly given those roots. The default of
#' [matern.operators()] and [spde.matern.operators()] is `"brasil"`, which is
#' right for `type = "covariance"`; a caller who never named a type should not
#' be refused because of it, so an untouched default becomes `"chebfunLB"`
#' here. A type the caller did name is passed through and refused if it is one
#' of the two.
#' @param type The value of `type_rational_approximation`.
#' @param user_set Did the caller name a type?
#' @return A single type name.
#' @noRd
wl2_operator_type <- function(type, user_set) {
  if (!user_set) {
    return("chebfunLB")
  }
  type[[1]]
}
