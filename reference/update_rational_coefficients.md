# Recompute the weighted-L2 coefficients for an estimated kappa

Two-stage use of the weighted-L2 rational coefficients. The coefficients
must not depend on `kappa` during estimation, since that would make the
likelihood non-smooth in `kappa`, so they are computed at set-up for a
conservative reference `kappa`. Once `kappa` has been estimated, this
function recomputes them once for a reference derived from the estimate,
for the final likelihood evaluation and for prediction.

## Usage

``` r
update_rational_coefficients(object, kappa_ref = NULL, safety = 3, ...)
```

## Arguments

- object:

  A model created by
  [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
  with `type_rational_approximation = "wl2"`.

- kappa_ref:

  The new reference value. If `NULL`, the `kappa` of the object divided
  by `safety` is used.

- safety:

  The safety factor applied to the `kappa` of the object when
  `kappa_ref` is not given.

- ...:

  Further arguments passed to
  [`update()`](https://rdrr.io/r/stats/update.html).

## Value

The updated model.

## Details

The reference must remain a lower bound for `kappa`, which is why the
default divides the estimate by `safety`. A reference above the true
`kappa` leaves part of the spectrum outside the fitted interval and the
error grows by one to two orders of magnitude, whereas a reference below
it degrades gracefully towards the mesh-free fit.

## See also

[`rational.coefficients.wl2()`](https://davidbolin.github.io/rSPDE/reference/rational.coefficients.wl2.md),
[`rspde.xmin()`](https://davidbolin.github.io/rSPDE/reference/rspde.xmin.md)

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 101)
op <- matern.operators(
  loc_mesh = x, nu = 0.4, range = 0.2, sigma = 1, d = 1, m = 2,
  parameterization = "matern", type = "operator",
  type_rational_approximation = "wl2"
)
op <- update_rational_coefficients(op)
```
