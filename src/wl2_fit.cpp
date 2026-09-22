/* Weighted-L2 rational coefficients: the inner fit in C++.
 *
 * Mirrors wl2_nnls(), wl2_resjac() and wl2_lm() of R/rational_wl2.R exactly.
 * The R side keeps the starts, the continuation in alpha and in m, the table
 * assembly and the cache; only the Levenberg-Marquardt loop over the poles,
 * with its non-negative least squares inner solve and Golub-Pereyra Jacobian,
 * lives here. That loop is where the time goes: the matrices are small
 * (n <= 900 rows, at most 8 columns) so in R it runs at a few per cent of the
 * machine's arithmetic throughput, dominated by interpreter and S4 dispatch
 * overhead rather than by flops.
 *
 * The standard library comes first: Rinternals.h defines length() as a macro,
 * which breaks <locale> and hence <vector> on libc++.
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

namespace {

const double kEps = 2.220446049250313e-16;

/* Least squares on the first p columns of A (n x p, column major, destroyed).
 * Rank deficient columns get coefficient zero, as .lm.fit(tol = 1e-14) does.
 * On return qr holds the factorisation and tau the reflectors, so the caller
 * can reuse them. */
bool ls_solve(const double *A, int n, int p, const double *b,
              std::vector<double> &coef, std::vector<double> &qr,
              std::vector<double> &tau, std::vector<double> &work) {
  qr.assign(A, A + (size_t)n * p);
  tau.assign(p, 0.0);
  int lwork = (int)work.size();
  int info = 0;
  F77_CALL(dgeqrf)(&n, &p, qr.data(), &n, tau.data(), work.data(), &lwork,
                   &info);
  if (info != 0) return false;

  std::vector<double> qtb(b, b + n);
  int one = 1;
  char side = 'L', trans = 'T';
  F77_CALL(dormqr)(&side, &trans, &n, &one, &p, qr.data(), &n, tau.data(),
                   qtb.data(), &n, work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) return false;

  /* back substitution, skipping negligible pivots */
  double rmax = 0.0;
  for (int j = 0; j < p; j++) {
    rmax = std::max(rmax, std::fabs(qr[(size_t)j * n + j]));
  }
  const double rtol = 1e-14 * (rmax > 0.0 ? rmax : 1.0);
  coef.assign(p, 0.0);
  for (int j = p - 1; j >= 0; j--) {
    const double d = qr[(size_t)j * n + j];
    if (std::fabs(d) <= rtol) {
      coef[j] = 0.0;
      continue;
    }
    double s = qtb[j];
    for (int k = j + 1; k < p; k++) s -= qr[(size_t)k * n + j] * coef[k];
    coef[j] = s / d;
  }
  return true;
}

/* Lawson-Hanson non-negative least squares on column-normalised columns. The
 * normalisation is essential: the columns span several orders of magnitude,
 * and the optimality test compares gradients across them. */
void nnls(const double *A, int n, int M, const double *f,
          std::vector<double> &x, std::vector<double> &scl,
          std::vector<double> &An, std::vector<double> &work) {
  scl.assign(M, 1.0);
  An.assign((size_t)n * M, 0.0);
  for (int j = 0; j < M; j++) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += A[(size_t)j * n + i] * A[(size_t)j * n + i];
    s = std::sqrt(s);
    if (!(s > 0.0)) s = 1.0;
    scl[j] = s;
    for (int i = 0; i < n; i++) An[(size_t)j * n + i] = A[(size_t)j * n + i] / s;
  }

  x.assign(M, 0.0);
  std::vector<char> P(M, 0);
  std::vector<double> w(M, 0.0), resid(f, f + n), s(M, 0.0);
  std::vector<double> sub((size_t)n * M), coef, qr, tau;

  for (int j = 0; j < M; j++) {
    double g = 0.0;
    for (int i = 0; i < n; i++) g += An[(size_t)j * n + i] * f[i];
    w[j] = g;
  }
  double wmax0 = 0.0;
  for (int j = 0; j < M; j++) wmax0 = std::max(wmax0, std::fabs(w[j]));
  const double tol = 10.0 * kEps * std::max(n, M) * wmax0;

  const int maxit = 10 * M;
  for (int it = 0; it < maxit; it++) {
    int best = -1;
    double bestw = tol;
    for (int j = 0; j < M; j++) {
      if (!P[j] && w[j] > bestw) {
        bestw = w[j];
        best = j;
      }
    }
    if (best < 0) break;
    P[best] = 1;

    for (int inner = 0; inner <= maxit; inner++) {
      int np = 0;
      std::vector<int> idx;
      for (int j = 0; j < M; j++) {
        if (P[j]) {
          for (int i = 0; i < n; i++) {
            sub[(size_t)np * n + i] = An[(size_t)j * n + i];
          }
          idx.push_back(j);
          np++;
        }
      }
      if (np == 0) break;
      if (!ls_solve(sub.data(), n, np, f, coef, qr, tau, work)) break;
      std::fill(s.begin(), s.end(), 0.0);
      for (int k = 0; k < np; k++) s[idx[k]] = coef[k];

      double smin = 0.0;
      bool first = true;
      for (int k = 0; k < np; k++) {
        if (first || s[idx[k]] < smin) {
          smin = s[idx[k]];
          first = false;
        }
      }
      if (smin >= 0.0) {
        x = s;
        break;
      }
      double step = 1.0;
      for (int k = 0; k < np; k++) {
        const int j = idx[k];
        if (s[j] < 0.0) {
          const double den = x[j] - s[j];
          if (den != 0.0) step = std::min(step, x[j] / den);
        }
      }
      double xmax = 0.0;
      for (int j = 0; j < M; j++) {
        x[j] += step * (s[j] - x[j]);
        xmax = std::max(xmax, std::fabs(x[j]));
      }
      const double xtol = 10.0 * kEps * xmax;
      bool any = false;
      for (int j = 0; j < M; j++) {
        if (P[j] && !(x[j] > xtol)) P[j] = 0;
        if (!P[j]) x[j] = 0.0;
        if (P[j]) any = true;
      }
      if (!any) break;
    }

    for (int i = 0; i < n; i++) {
      double v = 0.0;
      for (int j = 0; j < M; j++) v += An[(size_t)j * n + i] * x[j];
      resid[i] = f[i] - v;
    }
    for (int j = 0; j < M; j++) {
      double g = 0.0;
      for (int i = 0; i < n; i++) g += An[(size_t)j * n + i] * resid[i];
      w[j] = g;
    }
  }
  for (int j = 0; j < M; j++) x[j] /= scl[j];
}

struct Workspace {
  std::vector<double> U, g, b, A, c, r, J, scl, An, work, sub, coef, qr, tau, Q;
};

/* Residual, and optionally the Golub-Pereyra Jacobian. */
void resjac(const double *x, const double *sw, const double *f, int n,
            const double *theta, int ntheta, bool shifted, bool kterm,
            bool need_jac, Workspace &ws) {
  const int M = shifted ? ntheta - 1 : ntheta;
  /* One extra basis column, constant in x, when a constant term is allowed.
     It has no parameter of its own, so its Jacobian column is zero. */
  const int ncol = M + (kterm ? 1 : 0);
  ws.b.assign(M, 0.0);
  for (int j = 0; j < M; j++) {
    double t = theta[shifted ? j + 1 : j];
    if (t < -5.0) t = -5.0;
    if (t > 60.0) t = 60.0;
    ws.b[j] = std::exp(t);
  }
  ws.U.assign((size_t)n * M, 0.0);
  for (int j = 0; j < M; j++) {
    const double bm1 = ws.b[j] - 1.0;
    for (int i = 0; i < n; i++) {
      ws.U[(size_t)j * n + i] = x[i] / (1.0 + bm1 * x[i]);
    }
  }
  double b0 = 0.0;
  ws.g.assign(n, 1.0);
  if (shifted) {
    b0 = 1.0 / (1.0 + std::exp(-theta[0]));
    for (int i = 0; i < n; i++) ws.g[i] = x[i] / (1.0 + (b0 - 1.0) * x[i]);
  }
  ws.A.assign((size_t)n * ncol, 0.0);
  for (int j = 0; j < M; j++) {
    for (int i = 0; i < n; i++) {
      ws.A[(size_t)j * n + i] = ws.U[(size_t)j * n + i] * sw[i] *
                                (shifted ? ws.g[i] : 1.0);
    }
  }
  if (kterm) {
    for (int i = 0; i < n; i++) ws.A[(size_t)M * n + i] = sw[i];
  }

  nnls(ws.A.data(), n, ncol, f, ws.c, ws.scl, ws.An, ws.work);

  ws.r.assign(n, 0.0);
  for (int i = 0; i < n; i++) {
    double v = 0.0;
    for (int j = 0; j < ncol; j++) v += ws.A[(size_t)j * n + i] * ws.c[j];
    ws.r[i] = v - f[i];
  }
  if (!need_jac) return;

  std::vector<int> keep;
  for (int j = 0; j < ncol; j++) {
    if (ws.c[j] > 0.0) keep.push_back(j);
  }
  ws.J.assign((size_t)n * ntheta, 0.0);
  const int np = (int)keep.size();
  if (np == 0) return;

  /* QR of the passive columns, then form Q explicitly and reuse it */
  ws.sub.assign((size_t)n * np, 0.0);
  for (int k = 0; k < np; k++) {
    for (int i = 0; i < n; i++) {
      ws.sub[(size_t)k * n + i] = ws.A[(size_t)keep[k] * n + i];
    }
  }
  ws.qr = ws.sub;
  ws.tau.assign(np, 0.0);
  int lwork = (int)ws.work.size();
  int info = 0;
  F77_CALL(dgeqrf)(&n, &np, ws.qr.data(), &n, ws.tau.data(), ws.work.data(),
                   &lwork, &info);
  if (info != 0) return;
  std::vector<double> R((size_t)np * np, 0.0);
  for (int j = 0; j < np; j++) {
    for (int i = 0; i <= j; i++) R[(size_t)j * np + i] = ws.qr[(size_t)j * n + i];
  }
  ws.Q = ws.qr;
  F77_CALL(dorgqr)(&n, &np, &np, ws.Q.data(), &n, ws.tau.data(),
                   ws.work.data(), &lwork, &info);
  if (info != 0) return;

  std::vector<double> D((size_t)n * ncol), v(n), qtv(np), z(np), y(np);
  for (int k = 0; k < ntheta; k++) {
    std::fill(D.begin(), D.end(), 0.0);
    if (shifted && k == 0) {
      const double db0 = b0 * (1.0 - b0);
      for (int j = 0; j < M; j++) {
        for (int i = 0; i < n; i++) {
          D[(size_t)j * n + i] =
              -ws.g[i] * ws.g[i] * db0 * ws.U[(size_t)j * n + i] * sw[i];
        }
      }
    } else {
      const int j = shifted ? k - 1 : k;
      for (int i = 0; i < n; i++) {
        const double u = ws.U[(size_t)j * n + i];
        D[(size_t)j * n + i] =
            (shifted ? ws.g[i] : 1.0) * (-ws.b[j] * u * u) * sw[i];
      }
    }
    for (int i = 0; i < n; i++) {
      double s = 0.0;
      for (int j = 0; j < ncol; j++) s += D[(size_t)j * n + i] * ws.c[j];
      v[i] = s;
    }
    /* t1 = v - Q Q^T v */
    for (int k2 = 0; k2 < np; k2++) {
      double s = 0.0;
      for (int i = 0; i < n; i++) s += ws.Q[(size_t)k2 * n + i] * v[i];
      qtv[k2] = s;
    }
    /* z = D_P^T (-r) */
    for (int k2 = 0; k2 < np; k2++) {
      double s = 0.0;
      for (int i = 0; i < n; i++) {
        s += D[(size_t)keep[k2] * n + i] * (-ws.r[i]);
      }
      z[k2] = s;
    }
    /* y = R^{-T} z */
    for (int i = 0; i < np; i++) {
      double s = z[i];
      for (int k2 = 0; k2 < i; k2++) s -= R[(size_t)i * np + k2] * y[k2];
      const double d = R[(size_t)i * np + i];
      y[i] = (std::fabs(d) > 0.0) ? s / d : 0.0;
    }
    for (int i = 0; i < n; i++) {
      double t1 = v[i], t2 = 0.0;
      for (int k2 = 0; k2 < np; k2++) {
        t1 -= ws.Q[(size_t)k2 * n + i] * qtv[k2];
        t2 += ws.Q[(size_t)k2 * n + i] * y[k2];
      }
      ws.J[(size_t)k * n + i] = t1 + t2;
    }
  }
}

double sumsq(const std::vector<double> &v) {
  double s = 0.0;
  for (size_t i = 0; i < v.size(); i++) s += v[i] * v[i];
  return s;
}

} /* namespace */

extern "C" SEXP rspde_wl2_lm(SEXP x_, SEXP sw_, SEXP f_, SEXP theta0_,
                             SEXP shifted_, SEXP kterm_, SEXP ctrl_) {
  const int n = LENGTH(x_);
  const double *x = REAL(x_), *sw = REAL(sw_), *f = REAL(f_);
  const int ntheta = LENGTH(theta0_);
  const bool shifted = (LOGICAL(shifted_)[0] == TRUE);
  const bool kterm = (LOGICAL(kterm_)[0] == TRUE);
  const double ftol = REAL(ctrl_)[0], xtol = REAL(ctrl_)[1];
  const int maxfev = (int)REAL(ctrl_)[2], maxiter = (int)REAL(ctrl_)[3];

  std::vector<double> theta(REAL(theta0_), REAL(theta0_) + ntheta);
  Workspace ws;
  ws.work.assign(std::max(64 * (ntheta + 1), 256), 0.0);

  resjac(x, sw, f, n, theta.data(), ntheta, shifted, kterm, true, ws);
  double fval = sumsq(ws.r);
  int nfev = 1;
  double lambda = 1e-3;
  std::vector<double> J = ws.J;

  const int np = ntheta;
  std::vector<double> JtJ((size_t)np * np), Jtr(np), dg(np), Amat((size_t)np * np),
      step(np), theta_new(np);

  for (int iter = 0; iter < maxiter; iter++) {
    for (int a = 0; a < np; a++) {
      for (int b2 = 0; b2 < np; b2++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
          s += J[(size_t)a * n + i] * J[(size_t)b2 * n + i];
        }
        JtJ[(size_t)b2 * np + a] = s;
      }
      double s = 0.0;
      for (int i = 0; i < n; i++) s += J[(size_t)a * n + i] * ws.r[i];
      Jtr[a] = s;
    }
    double dgmax = 0.0;
    for (int a = 0; a < np; a++) dgmax = std::max(dgmax, JtJ[(size_t)a * np + a]);
    const double dgfloor = 1e-10 * std::max(dgmax, 1e-300);
    for (int a = 0; a < np; a++) {
      dg[a] = std::max(JtJ[(size_t)a * np + a], dgfloor);
    }
    if (nfev > maxfev) break;

    bool accepted = false;
    while (lambda < 1e14) {
      Amat = JtJ;
      for (int a = 0; a < np; a++) Amat[(size_t)a * np + a] += lambda * dg[a];
      for (int a = 0; a < np; a++) step[a] = -Jtr[a];
      int info = 0, one = 1;
      char uplo = 'U';
      F77_CALL(dposv)(&uplo, &np, &one, Amat.data(), &np, step.data(), &np,
                      &info FCONE);
      bool ok = (info == 0);
      if (ok) {
        for (int a = 0; a < np; a++) {
          if (!R_FINITE(step[a])) ok = false;
        }
      }
      if (!ok) {
        lambda *= 10.0;
        continue;
      }
      for (int a = 0; a < np; a++) theta_new[a] = theta[a] + step[a];
      resjac(x, sw, f, n, theta_new.data(), ntheta, shifted, kterm, true, ws);
      const double fnew = sumsq(ws.r);
      nfev++;
      if (R_FINITE(fnew) && fnew < fval) {
        const double df = fval - fnew;
        double dx = 0.0, tn = 0.0;
        for (int a = 0; a < np; a++) {
          dx += step[a] * step[a];
          tn += theta_new[a] * theta_new[a];
        }
        dx = std::sqrt(dx);
        tn = std::sqrt(tn);
        theta = theta_new;
        J = ws.J;
        fval = fnew;
        lambda = std::max(lambda / 10.0, 1e-12);
        accepted = true;
        if (df <= ftol * fval || dx <= xtol * (tn + xtol)) {
          iter = maxiter; /* converged */
        }
        break;
      }
      lambda *= 10.0;
      if (nfev > maxfev) break;
    }
    if (!accepted || nfev > maxfev) break;
  }

  /* final state at the returned theta */
  resjac(x, sw, f, n, theta.data(), ntheta, shifted, kterm, false, ws);

  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP th = PROTECT(allocVector(REALSXP, ntheta));
  for (int a = 0; a < ntheta; a++) REAL(th)[a] = theta[a];
  SEXP cost = PROTECT(ScalarReal(sumsq(ws.r) / 2.0));
  SEXP cc = PROTECT(allocVector(REALSXP, (int)ws.c.size()));
  for (size_t j = 0; j < ws.c.size(); j++) REAL(cc)[j] = ws.c[j];
  SET_VECTOR_ELT(out, 0, th);
  SET_VECTOR_ELT(out, 1, cost);
  SET_VECTOR_ELT(out, 2, cc);
  SEXP nms = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(nms, 0, mkChar("theta"));
  SET_STRING_ELT(nms, 1, mkChar("cost"));
  SET_STRING_ELT(nms, 2, mkChar("c"));
  setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(5);
  return out;
}

/* ------------------------------------------------------------------ */

extern "C" {

static const R_CallMethodDef CallEntries[] = {
    {"rspde_wl2_lm", (DL_FUNC)&rspde_wl2_lm, 7},
    {NULL, NULL, 0}};

void R_init_rSPDE(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

} /* extern "C" */
