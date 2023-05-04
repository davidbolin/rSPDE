# Description #

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/rSPDE)](https://cran.r-project.org/package=rSPDE)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)
[![R-CMD-check](https://github.com/davidbolin/rSPDE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidbolin/rSPDE/actions/workflows/R-CMD-check.yaml)

`rSPDE` is an R package used for computing rational approximations of fractional SPDEs. These rational approximations can be used for computationally efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented. The package also contains interfaces to [R-INLA][ref4] and [inlabru][ref8].

For illustration purposes, the package contains a simple FEM implementation for models on `R`. See the 
[Getting Started to the rSPDE package][ref2] vignette for an introduction to the package. The [Rational approximation with the rSPDE package][ref6] and [Operator-based rational approximation ][ref5] vignettes provide
introductions to how to create and fit rSPDE models. For an introduction to the R-INLA implementation
of the rSPDE models see the [R-INLA implementation of the rational SPDE approach][ref3]. The [`rSPDE` documentation][ref7] contains descriptions and examples of the functions in the `rSPDE` package.

# Shiny app

The paper [Covariance-based rational approximations of fractional SPDEs for computationally efficient Bayesian inference][ref9] contains a study on the accuracy of the various rational approximations that are implemented in the `rSPDE` package. These are summarized in a shiny app that can be run directly from GitHub by using the following command (you may need to install some of the dependencies):

```r
# install.packages("shiny", "shinythemes", "plotly")
# install.packages("tidyr", "dplyr")
shiny::runGitHub("davidbolin/rSPDE", subdir="shiny_app")
```

# References #
Z. Xiong, A. Simas, D. Bolin (2022) [Covariance-based rational approximations of fractional SPDEs for computationally efficient Bayesian inference][ref9]. 	ArXiv:2209.04670

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
where the variables `rtools` and `gcc` need to be changed if `Rtools` is not installed directly on `C:`,
and `gcc`'s version might need to be changed depending on the version of `Rtools`.

Finally, if you want to build the `rSPDE` package from source on Mac or Linux, please
click [here](https://davidbolin.github.io/rSPDE//articles/build_source.html).

# Repository branch workflows #
The package version format for released versions is `major.minor.bugfix`. All regular development should be performed on the `devel` branch or in a feature branch, managed with `git flow feature`. Ideally, all the changes should be made on the `devel` branch. The `devel` version of the package should contain unit tests and examples for all important functions. Several functions may depend on `INLA`. Examples and tests for such functions might create problems when submitting to CRAN. To solve this problem, we created some Github Actions scripts that get the examples and tests depending on `INLA` on the `devel` branch and adapt to versions that will not fail on CRAN. Therefore, the best way to handle these situations is to avoid as much as possible to do any push to the `stable` branch. The idea is to update the `stable` branch by merges following the workflow that will be described below. 
The examples that depend on `INLA` should have the following structure:

```
#' \donttest{ #devel version
#' library(INLA)
#' 
#' # The contents of the example
#'
#' #devel.tag
#' }
```

The tests that depend on `INLA` should have the following structure:

```
test_that("Description of the test", {
  testthat::skip_on_cran()
  inla_installed <- "INLA" %in% rownames(installed.packages())
  if(!inla_installed){
    testthat::skip("INLA not installed")
  }
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  
  # The contents of the test
  
  INLA::inla.setOption(num.threads = old_threads)
})
```

On the `devel` branch, the vestion number is `major.minor.bugfix.9000`, where the first three components reflect the latest released version with changes present in the `default` branch. Bugfixes should be applied via the `git flow bugfix` and `git flow hotfix` methods, as indicated below. For `git flow` configuration, use `master` as the stable master branch, `devel` as the develop branch, and `v` as the version tag prefix. Hotfixes directly `stable` should be avoided whenever possible to minimize conflicts on merges. See [the `git flow` tutorial](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more information.

For non `master` and `devel` branches that collaborators need access to (e.g. release branches, feature branches, etc, use the `git flow publish` mechanism).


  * Prepare a new stable release with CRAN submission:
```
git flow release start major.(minor+1).0

# In R, the following updates the version number in DESCRIPTION and NEWS:
usethis::use_version("minor") 
## At this point, see the CRAN submission section below.
git flow release finish 'VERSION'
# In the stable merge, accept all incoming changes.
# Push the changes and do adjustments if needed.
# After handling the stable branch, merge back with devel.
# In R, the following updates the dev version number in DESCRIPTION and NEWS:
usethis::use_dev_version() 
```
  * Do a hotfix (branch from stable branch; use bugfix for release branch bugfixes):
```
git flow hotfix start hotfix_branch_name
## Do the bugfix, update the verison number major.minor.(bugfix+1), and commit
## Optionally, do CRAN submission
git flow hotfix finish hotfix_branch_name
## Resolve merge conflicts, if any
```
  * CRAN submission
```
## Perform CRAN checks on the stable branch
## If unsuccessful then do bugfixes with increasing bugfix version, until ok
## Submit to CRAN
## If not accepted then do more bugfixes and repeat
```



[ref]: https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1665537  "The rational SPDE approach for Gaussian random fields with general smoothness"
[ref2]: https://davidbolin.github.io/rSPDE//articles/rSPDE.html "Getting Started to the rSPDE package"
[ref3]: https://davidbolin.github.io/rSPDE//articles/rspde_inla.html "INLA Vignette"
[ref4]: https://r-inla.org "INLA homepage"
[ref5]: https://davidbolin.github.io/rSPDE//articles/rspde_base.html
[ref6]: https://davidbolin.github.io/rSPDE//articles/rspde_cov.html
[ref7]: https://davidbolin.github.io/rSPDE/reference/index.html "`rSPDE` documentation."
[ref8]: https://sites.google.com/inlabru.org/inlabru "inlabru homepage"
[ref9]: https://arxiv.org/abs/2209.04670 "Covariance-based rational approximations of fractional SPDEs for computationally efficient Bayesian inference"
