# Lower end of the spectral interval of a rational approximation

Computes the value \\x\_{\min} = 1/(1 +
\mu\_{\max}/\kappa\_{\mathrm{lo}}^2)\\ used by the weighted-L2 rational
coefficients, where \\\mu\_{\max}\\ is the largest generalised
eigenvalue of \\(G, C)\\ and \\\kappa\_{\mathrm{lo}}\\ is a lower bound
for the values that \\\kappa\\ can take.

## Usage

``` r
rspde.xmin(
  C = NULL,
  G = NULL,
  mesh = NULL,
  kappa_ref = NULL,
  nu = NULL,
  diameter = NULL,
  loc_mesh = NULL,
  eigenvalue = c("bound", "exact")
)
```

## Arguments

- C:

  The mass matrix of the finite element discretisation.

- G:

  The stiffness matrix of the finite element discretisation.

- mesh:

  An optional mesh; `C` and `G` are computed from it if they are not
  given, and it is used for the diameter of the domain.

- kappa_ref:

  The reference value \\\kappa\_{\mathrm{lo}}\\. If `NULL`, it is taken
  to be `sqrt(8 * nu) / diameter`.

- nu:

  The smoothness parameter, used for the default `kappa_ref`.

- diameter:

  The diameter of the domain, used for the default `kappa_ref`. Computed
  from `mesh` if not given.

- loc_mesh:

  Mesh locations, an alternative to `mesh` for the diameter.

- eigenvalue:

  Either `"bound"` (the default), which uses the Gershgorin bound
  `max(rowSums(abs(G)) / diag(C))` for \\\mu\_{\max}\\, or `"exact"`,
  which computes the eigenvalue by Lanczos iteration. The bound is an
  over-estimate of \\\mu\_{\max}\\, and therefore gives a conservative
  (too small) \\x\_{\min}\\, which is the safe direction.

## Value

The value of \\x\_{\min}\\.

## Details

The reference \\\kappa\_{\mathrm{lo}}\\ must be a *lower* bound: a
reference above the true \\\kappa\\ leaves part of the spectrum outside
the fitted interval, and the error then grows by one to two orders of
magnitude. A reference below the true \\\kappa\\ is safe; the
approximation then degrades gracefully towards the mesh-free fit. The
default, \\\kappa\_{\mathrm{lo}} = \sqrt{8\nu}/\mathrm{diam}\\,
corresponds to a range equal to the diameter of the domain.

## See also

[`rational.coefficients.wl2()`](https://davidbolin.github.io/rSPDE/reference/rational.coefficients.wl2.md)

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 201)
fem <- rSPDE.fem1d(x)
# reference kappa corresponding to a range equal to the domain diameter
rspde.xmin(C = fem$C, G = fem$G, loc_mesh = x, nu = 0.5)
#> [1] 2.499938e-05
# a user-supplied lower bound for kappa
rspde.xmin(C = fem$C, G = fem$G, kappa_ref = 10)
#> [1] 0.0006246096
```
