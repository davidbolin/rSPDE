# Warnings free loading of add-on packages

Turn off all warnings for require(), to allow clean completion of
examples that require unavailable Suggested packages.

## Usage

``` r
require.nowarnings(package, lib.loc = NULL, character.only = FALSE)
```

## Arguments

- package:

  The name of a package, given as a character string.

- lib.loc:

  a character vector describing the location of R library trees to
  search through, or `NULL`. The default value of `NULL` corresponds to
  all libraries currently known to
  [`.libPaths()`](https://rdrr.io/r/base/libPaths.html). Non-existent
  library trees are silently ignored.

- character.only:

  a logical indicating whether `package` can be assumed to be a
  character string.

## Value

`require.nowarnings` returns (invisibly) `TRUE` if it succeeds,
otherwise `FALSE`

## Details

[`require(package)`](https://rdrr.io/r/base/library.html) acts the same
as
[`require(package, quietly = TRUE)`](https://rdrr.io/r/base/library.html)
but with warnings turned off. In particular, no warning or error is
given if the package is unavailable. Most cases should use
[`requireNamespace(package, quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html)
instead, which doesn't produce warnings.

## See also

[`require()`](https://rdrr.io/r/base/library.html)

## Examples

``` r
## This should produce no output:
if (require.nowarnings(nonexistent)) {
  message("Package loaded successfully")
}
```
