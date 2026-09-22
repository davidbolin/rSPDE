# Weighted-L2 rational coefficients

Computes the coefficients of a rational approximation of \\x^{\alpha}\\
(covariance type) or \\x^{\alpha/2}\\ (operator type) on the interval
\\\[x\_{\min}, 1\]\\, minimising the \\L_2\\ error weighted by the Weyl
density \\w_d(x) = x^{-1-d/2}(1-x)^{d/2-1}\\.

## Usage

``` r
rational.coefficients.wl2(
  alpha,
  d,
  m,
  x_min = NULL,
  type = c("covariance", "operator"),
  n_starts = 8,
  n_grid = 300,
  start = NULL,
  seed = 1L,
  continuation = TRUE,
  s = 0,
  k_term = FALSE
)
```

## Arguments

- alpha:

  The exponent, \\\alpha = \nu + d/2\\. For the covariance type,
  \\\lfloor \alpha\rfloor\\ must be 0, 1 or 2; for the operator type,
  \\\alpha\\ must be smaller than 2.

- d:

  The dimension of the domain.

- m:

  The order of the rational approximation.

- x_min:

  Lower end of the spectral interval, or `NULL` for the mesh-free fit on
  \\(0,1\]\\. See
  [`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md).

- type:

  Either `"covariance"` or `"operator"`.

- n_starts:

  Number of starting values for the outer optimisation. Eight starts are
  needed for `m = 4`; for `m` at most 3 a single start suffices.

- n_grid:

  Number of quadrature nodes.

- start:

  Optional starting value for the outer parameters, for instance the
  `theta` element returned by a fit for a neighbouring value of `alpha`.

- seed:

  Seed of the local random number stream used for the random starting
  values. The user's random number stream is not affected.

- continuation:

  Should the fit use continuation in `m`, solving for `1, ..., m` and
  starting each order from the poles of the previous one with one pole
  inserted? This is needed beyond `m = 4`, where a plain multistart
  tends to settle at the optimum of the next smaller order. Ignored when
  `start` is supplied.

- s:

  Exponent of an extra factor \\x^s\\ in the weight, on top of the Weyl
  density. `s = 0`, the default, is the weight for which the objective
  is the \\L_2\\ error of the covariance. Larger `s` de-emphasises the
  small eigenvalues, and corresponds to measuring the error of the
  solution operator against data in \\H^s\\ rather than in \\L_2\\.

- k_term:

  Should a constant term be included, as in the tabulated coefficients?
  It is fitted, non-negative, and returned as `k`. The default `FALSE`
  is what the covariance-based models need, since a constant term makes
  the trace infinite. It only makes sense together with a finite
  `x_min`: on \\(0,1\]\\ with `s` below \\d/2\\ the objective is
  infinite for any non-zero constant, and the fit drives it to zero.

## Value

A list with elements

- r:

  The residues, all non-negative.

- p:

  The poles, all smaller than one.

- p0:

  The shift of the integer factor, or `NULL` if \\\lfloor\alpha\rfloor =
  0\\.

- q:

  The power of the integer factor, \\\lfloor\alpha\rfloor\\ for the
  covariance type and 0 for the operator type.

- k:

  The constant term: 0 unless `k_term` is `TRUE`.

- rel_err:

  The relative weighted \\L_2\\ error of the fit.

- x_min:

  The lower end of the interval that was used.

- theta:

  The internal parameters, for warm starts.

- kind:

  `"plain"`, `"shifted"` or `"shifted2"`, for `q` equal to 0, 1 and 2.

## Details

Here \\x = 1/\lambda\\, where \\\lambda\\ are the eigenvalues of \\L_h =
\kappa^{-2}(\kappa^2 C + G)\\ relative to \\C\\, so that \\\lambda \ge
1\\ and \\x \in (0, 1\]\\. The value \\x\_{\min}\\ corresponds to the
largest eigenvalue, i.e. to the resolution of the mesh; `x_min = NULL`
gives the mesh-free fit on \\(0,1\]\\.

The classes are parameterised so that the resulting model is guaranteed
to be valid (all residues are non-negative and all poles are smaller
than one). For the covariance type the class carries an integer factor
of order \\q = \lfloor\alpha\rfloor\\, with one shift \\p_0\\ shared by
its factors:

- covariance:

  \\r(x) = (x/(1-p_0 x))^q\sum\_{j=1}^m r_j x/(1-p_j x)\\, for \\q =
  0\\, 1 and 2, the factor being absent when \\q = 0\\.

- operator:

  \\r(x) = \sum\_{i=1}^{m+1} r_i x/(1-p_i x)\\, approximating
  \\x^{\alpha/2}\\.

In contrast to the tabulated `"brasil"`, `"chebfun"` and `"chebfunLB"`
coefficients, there is no constant term, so that the covariance-based
models built from these coefficients have \\m\\ instead of \\m+1\\
blocks.

The fit is a variable projection: the poles are found by
Levenberg-Marquardt, with the residues eliminated by non-negative least
squares at every step. This inner loop is compiled (`src/wl2_fit.cpp`)
and is a few times faster than the equivalent R code, which is kept as a
reference implementation and is used instead when the option
`rSPDE.wl2.use.cpp` is set to `FALSE`.

## See also

[`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
cf <- rational.coefficients.wl2(alpha = 0.75, d = 1, m = 2)
cf$rel_err
#> [1] 0.01426543
```
