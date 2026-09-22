# Where generated coefficient tables are kept between sessions

The mesh-free weighted-L2 tables are shipped with the package, but a
table for a particular spectral interval (see
[`rspde.wl2.table()`](https://davidbolin.github.io/rSPDE/reference/rspde.wl2.table.md)
and the `kappa_ref` argument of
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md))
has to be computed, which takes of the order of a second. Within a
session such tables are cached in memory. Enabling this cache stores
them on disk as well, so that a second session, or a script run again,
finds them instead of refitting.

The cache is off by default: a package should not write outside the
session temporary directory unless asked to. Calling `rspde.cache(TRUE)`
is that request, and creates the directory. The environment variable
`RSPDE_CACHE_DIR` sets it for non-interactive use, for instance from
`.Renviron`.

Cached files are organised by package version, so upgrading the package
does not read tables produced by an older fit. Within a development
version, clear the cache after changing anything that affects the
coefficients.

## Usage

``` r
rspde.cache(dir, clear = FALSE)
```

## Arguments

- dir:

  `TRUE` to enable the cache at the standard location for this package,
  `tools::R_user_dir("rSPDE", "cache")`; `FALSE` to disable it; a path
  to use that directory instead; or missing to query the current
  setting.

- clear:

  Delete the tables that are currently stored?

## Value

The cache directory, or `NULL` if the cache is off, invisibly when
called for its effect.

## See also

[`rspde.wl2.table()`](https://davidbolin.github.io/rSPDE/reference/rspde.wl2.table.md),
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)

## Examples

``` r
rspde.cache()
#> NULL
if (FALSE) { # \dontrun{
rspde.cache(TRUE) # store generated tables between sessions
rspde.cache(clear = TRUE) # throw away what is stored
rspde.cache(FALSE)
} # }
```
