# Transform Spacetime SPDE Model Parameters to Original Scale

This function takes a vector of transformed parameters and applies the
appropriate transformations to return them in the original scale for use
in spacetime SPDE models.

## Usage

``` r
transform_parameters_spacetime(theta, st_model)
```

## Arguments

- theta:

  A numeric vector containing the transformed parameters in this order:

  lkappa

  :   The logarithmic representation of kappa.

  lsigma

  :   The logarithmic representation of sigma.

  lgamma

  :   The logarithmic representation of gamma.

  logit_rho (optional)

  :   The logit-transformed representation of rho, if drift = 1.

  logit_rho2 (optional)

  :   The logit-transformed representation of rho2, if drift = 1 and d =
      2.

- st_model:

  A list containing the spacetime model parameters:

  d

  :   The dimension (e.g., 1 or 2).

  bound

  :   The bound for rho and rho2.

  is_bounded

  :   A logical value indicating if rho and rho2 are bounded.

  drift

  :   A logical value indicating if drift is included in the model.

## Value

A named list with the parameters in the original scale:

- kappa:

  The original scale for kappa (exponential of lkappa).

- sigma:

  The original scale for sigma (exponential of lsigma).

- gamma:

  The original scale for gamma (exponential of lgamma).

- rho (optional):

  The original scale for rho.

- rho2 (optional):

  The original scale for rho2, if d = 2.
