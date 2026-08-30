# Load INLA safely for examples and tests

Checks that a sufficiently recent INLA package and matching executable
are available. In non-interactive use it also limits INLA to one thread
for repeatability. An optional Cgeneric symbol can be required when an
rSPDE feature depends on a model that may not yet be included in every
INLA build.

## Usage

``` r
rspde_safe_inla(
  multicore = NULL,
  quietly = FALSE,
  minimum_version = "24.12.01",
  required_symbol = NULL
)

local_rspde_safe_inla(
  multicore = FALSE,
  quietly = TRUE,
  required_symbol = NULL,
  envir = parent.frame()
)
```

## Arguments

- multicore:

  Logical; if `TRUE`, leave INLA's `num.threads` option unchanged. If
  `FALSE`, set it to `"1:1:1"`. The default allows multicore only in
  interactive sessions outside testthat.

- quietly:

  Logical; suppress diagnostic messages when `TRUE`.

- minimum_version:

  Minimum acceptable INLA version. This should match the requirement in
  `DESCRIPTION`.

- required_symbol:

  Optional name of a required INLA Cgeneric symbol.

- envir:

  Environment in which test cleanup handlers are registered.

## Value

`TRUE` when INLA is available and safe to use; otherwise `FALSE`.

`local_rspde_safe_inla()` is called for its testthat skip and local
option-setting side effects.

## Examples

``` r
if (rspde_safe_inla()) {
  # Run INLA-dependent calculations.
}
#> NULL
```
