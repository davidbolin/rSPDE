# Description #

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/rSPDE)](https://cran.r-project.org/package=rSPDE)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)](https://cranlogs.r-pkg.org/badges/grand-total/rSPDE)

rSPDE is an R package used for computing rational approximations of fractional SPDEs These rational approximations can be used for computatially efficient statistical inference.

Basic statistical operations such as likelihood evaluations and kriging predictions using the fractional approximations are also implemented.

For illustration purposes, the package contains a simple FEM implementation for models on R. See the 
[Vignette][ref2] for an introduction to the package. 

# Reference #
D. Bolin and K. Kichner, [The rational SPDE approach for Gaussian random fields with general smoothness][ref]. Journal of Computational and Graphical Statistics.

# Installation instructions #
The latest CRAN release of the package can be installed directly from CRAN with `install.packages("rSPDE")`.
The latest stable version (which is sometimes slightly more recent than the CRAN version), can be installed by using the command
```r
remotes::install_bitbucket("davidbolin/rSPDE", ref = "master")
```
in R. The development version can be installed using the command
```r
remotes::install_bitbucket("davidbolin/rSPDE", ref = "devel")
```


# Repository branch workflows #
The package version format for released versions is `major.minor.bugfix`. All regular development should be performed on the `devel` branch or in a feature branch, managed with `git flow feature`. On the `devel` branch, the vestion number is `major.minor.bugfix.9000`, where the first three components reflect the latest released version with changes present in the `default` branch. Bugfixes should be applied via the `git flow bugfix` and `git flow hotfix` methods, as indicated below. For `git flow` configuration, use `master` as the stable master branch, `devel` as the develop branch, and `v` as the version tag prefix. See [the `git flow` tutorial](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more information.

For non `master` and `devel` branches that collaborators need access to (e.g. release branches, feature branches, etc, use the `git flow publish` mechanism).


* Prepare a new stable release with CRAN submission:
```
git flow release start major.(minor+1).0
## Update the DESCRIPTION version number as major.(minor+1).0
## Update the version in NEWS.md
## Commit the changes
## At this point, see the CRAN submission section below.
git flow release finish 'VERSION'
## Resolve/update the DESCRIPTION and NEWS.md version number conflict
## in favour of the released version, with extra .9000, e.g. with
## the help of  git mergetool
## Add a new version section in NEWS.md
## Commit the merge
```

* Do a hotfix (branch from stable master; use bugfix for release branch bugfixes):
```
git flow hotfix start hotfix_branch_name
## Do the bugfix, update the verison number major.minor.(bugfix+1), and commit
## Optionally, do CRAN submission
git flow hotfix finish hotfix_branch_name
## Resolve merge conflicts (hopefully mostly due to version numbers)
```

* CRAN submission
```
## Perform CRAN checks (usually on the release branch version)
## If unsuccessful then do bugfixes with increasing bugfix version, until ok
## Submit to CRAN
## If not accepted then do more bugfixes and repeat
```


[ref]: https://www.tandfonline.com/doi/full/10.1080/10618600.2019.1665537  "The rational SPDE approach for Gaussian random fields with general smoothness"
[ref2]: https://cran.r-project.org/web/packages/rSPDE/vignettes/rspde.html "Vignette"