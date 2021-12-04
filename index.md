# `rSPDE` Package #

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/rSPDE)](https://cran.r-project.org/package=rSPDE)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)

`rSPDE` is an R package used for computing rational approximations of fractional SPDEs These rational approximations can be used for computatially efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented. The package also contains an interface to [R-INLA][ref4].

For illustration purposes, the package contains a simple FEM implementation for models on R. See the 
[Getting Started to the rSPDE package][ref2] vignette for an introduction to the package. The [Rational approximation][ref5] and [Covariance-based rational approximation][ref6] vignettes provide
introductions to how to create and fit rSPDE models. For an introduction to the R-INLA implementation
of the rSPDE models see the [INLA Vignette][ref3]. The [`rSPDE` documentation][ref7] contains description and examples of the functions in the `rSPDE` package.

# Reference #
D. Bolin and K. Kirchner (2020) [The rational SPDE approach for Gaussian random fields with general smoothness][ref]. Journal of Computational and Graphical Statistics, 29:2, 274-285.

# Installation instructions #
The latest CRAN release of the package can be installed directly from CRAN with `install.packages("rSPDE")`.
The latest stable version (which is sometimes slightly more recent than the CRAN version), can be installed by using the command
```r
remotes::install_github("davidbolin/rspde", ref = "stable")
```
in R. The development version can be installed using the command
```r
remotes::install_github("davidbolin/rspde", ref = "devel")
```

If you want to install the package using the `remotes::install_github`-method on Windows, you first need to install `Rtools` and add the paths to `Rtools` and `gcc` to the Windows `PATH` environment variable. This can be done for the current R session only using the commands
```r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH = paste(c(gcc, rtools, Sys.getenv("PATH")), collapse = ";"))
```
where the variables `rtools` and `gcc` need to be changed if `Rtool`s is not installed directly on `C:`.

# Upcoming features

- Implementation of covariance-based non-stationary models.
- Implementation of the covariance-based rational approximation on R-STAN interface.
- Implementation of covariance-based method to more general SPDE models.
- Implementation of PC-priors for R-INLA `rSPDE` models.

[ref]: https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1665537  "The rational SPDE approach for Gaussian random fields with general smoothness"
[ref2]: https://davidbolin.github.io/rSPDE/articles/rSPDE.html "Getting Started to the rSPDE package"
[ref3]: https://davidbolin.github.io/rSPDE/articles/inla_rspde "INLA Vignette"
[ref4]: https://r-inla.org "INLA homepage"
[ref5]: https://davidbolin.github.io/rSPDE/articles/rspde_base.html
[ref6]: https://davidbolin.github.io/rSPDE/articles/rspde_cov.html
[ref7]: https://davidbolin.github.io/rSPDE/reference/index.html "`rSPDE` documentation."
