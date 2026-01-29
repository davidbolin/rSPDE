# The 2d folded Matern covariance function

`folded.matern.covariance.2d` evaluates the 2d folded Matern covariance
function over an interval \\\[0,L\]\times \[0,L\]\\.

## Usage

``` r
folded.matern.covariance.2d(
  h,
  m,
  kappa,
  nu,
  sigma,
  L = 1,
  N = 10,
  boundary = c("neumann", "dirichlet", "periodic", "R2")
)
```

## Arguments

- h, m:

  Vectors with two coordinates.

- kappa:

  Range parameter.

- nu:

  Shape parameter.

- sigma:

  Standard deviation.

- L:

  The upper bound of the square \\\[0,L\]\times \[0,L\]\\. By default,
  `L=1`.

- N:

  The truncation parameter.

- boundary:

  The boundary condition. The possible conditions are `"neumann"`
  (default), `"dirichlet"`, `"periodic"` or `"R2"`.

## Value

The correspoding covariance.

## Details

`folded.matern.covariance.2d` evaluates the 1d folded Matern covariance
function over an interval \\\[0,L\]\times \[0,L\]\\ under different
boundary conditions. For periodic boundary conditions
\$\$C\_{\mathcal{P}}((h_1,h_2),(m_1,m_2)) = \sum\_{k_2=-\infty}^\infty
\sum\_{k_1=-\infty}^{\infty} (C(\\(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\\),\$\$
for Neumann boundary conditions
\$\$C\_{\mathcal{N}}((h_1,h_2),(m_1,m_2)) = \sum\_{k_2=-\infty}^\infty
\sum\_{k_1=-\infty}^{\infty}
(C(\\(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\\)+C(\\(h_1-m_1+2k_1L,
h_2+m_2+2k_2L)\\)+C(\\(h_1+m_1+2k_1L,h_2-m_2+2k_2L)\\)+
C(\\(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\\)),\$\$ and for Dirichlet boundary
conditions: \$\$C\_{\mathcal{D}}((h_1,h_2),(m_1,m_2)) =
\sum\_{k_2=-\infty}^\infty \sum\_{k_1=-\infty}^{\infty}
(C(\\(h_1-m_1+2k_1L,h_2-m_2+2k_2L)\\)-
C(\\(h_1-m_1+2k_1L,h_2+m_2+2k_2L)\\)-C(\\(h_1+m_1+2k_1L,
h_2-m_2+2k_2L)\\)+C(\\(h_1+m_1+2k_1L,h_2+m_2+2k_2L)\\)),\$\$ where
\\C(\cdot)\\ is the Matern covariance function: \$\$C(h) =
\frac{\sigma^2}{2^{\nu-1}\Gamma(\nu)} (\kappa h)^\nu K\_\nu(\kappa
h).\$\$

We consider the truncation for \\k_1,k_2\\ from \\-N\\ to \\N\\.

## Examples

``` r
h <- c(0.5, 0.5)
m <- c(0.5, 0.5)
folded.matern.covariance.2d(h, m, kappa = 10, nu = 1 / 5, sigma = 1)
#> [1] 1.000043
```
