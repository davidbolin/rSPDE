# Operations with the Pr and Pl operators

Functions for multiplying and solving with the \\P_r\\ and \\P_l\\
operators as well as the latent precision matrix \\Q = P_l C^{-1}P_l\\
and covariance matrix \\\Sigma = P_r Q^{-1} P_r^T\\. These operations
are done without first assembling \\P_r\\, \\P_l\\ in order to avoid
numerical problems caused by ill-conditioned matrices.

## Usage

``` r
Pr.mult(obj, v, transpose = FALSE)

Pr.solve(obj, v, transpose = FALSE)

Pl.mult(obj, v, transpose = FALSE)

Pl.solve(obj, v, transpose = FALSE)

Q.mult(obj, v)

Q.solve(obj, v)

Qsqrt.mult(obj, v, transpose = FALSE)

Qsqrt.solve(obj, v, transpose = FALSE)

Sigma.mult(obj, v)

Sigma.solve(obj, v)
```

## Arguments

- obj:

  rSPDE object

- v:

  vector to apply the operation to

- transpose:

  set to TRUE if the operation should be performed with the transposed
  object

## Value

A vector with the values of the operation

## Details

`Pl.mult`, `Pr.mult`, and `Q.mult` multiplies the vector with the
respective object. Changing `mult` to `solve` in the function names
multiplies the vector with the inverse of the object. `Qsqrt.mult` and
`Qsqrt.solve` performs the operations with the square-root type object
\\Q_r = C^{-1/2}P_l\\ defined so that \\Q = Q_r^T Q_r\\.
