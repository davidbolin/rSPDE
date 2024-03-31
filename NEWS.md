# rSPDE (development version)

* Created the `group_predict` function, to obtain predictions on a testing set based on observations on a training set.

# rSPDE 2.3.3
* Bugfix on rspde_lme when fitting with fixed smoothness.
* Added a 2d fem interface.
* Moved from using INLA's mesh functions to fmesher's mesh functions.
* Removing rgdal from suggests.
* The `data` argument in `predict.rspde_lme` has been changed to `newdata`.
* Adding `covariance_mesh` and `cov_function_mesh` methods as functions in the list returned by objects obtained from `matern.operators()` and `spde.matern.operators()`.
* Updated the internal structure to match the updates from the `MetricGraph` package.
* Updated the `cross_validation` function to match the updates in `inlabru`.
* Added `glance` and `augment` methods for `rspde_lme` objects.

# rSPDE 2.3.2
* Small improvement on speed for rspde_lme.
* Bugfix on Q for small values of nu in dimension 1.
* Adding parameterization option for rspde.result.
* Bugfix on which_repl in rspde_lme.
* Addressing issues related to the new version of the Matrix package.

# rSPDE 2.3.1
* Adding references in DESCRIPTION.
* Changing link to eigen library.

# rSPDE 2.3.0
* Fixed a bug on rSPDE.construct.matern.loglike when the parameterization is "matern".
* Created the rspde_lme() interface, with corresponding standard methods(predict, summary, etc).
* Updated the vignettes to use the rspde_lme() interface instead of the likelihood function factory.
* Replaced chol by Cholesky when using it to compute determinants or to solve systems.

# rSPDE 2.2.0
* Adding a new parameterization (variance and a range-like parameter)
* Posterior sampling on the predict method.
* Added the `cross_validation` function which has several scoring rules implemented (MSE, CRPS, SCRPS, DSS) based on our `inlabru` implementation of the rational SPDE approach.

# rSPDE 2.1.0
* Expanded the parameterization options on matern.operators and spde.matern.operators, along with their associated functions.
* Implementation of the precision method for inla_rspde objects.
* Implementation of the covariance-based spde.matern.operators function and its associated functions.
* Adjusts on the compatibility with the forthcoming MetricGraph package.

# rSPDE 2.0.0
* Added cgeneric versions of the nonstationary models
* Added support for metric graphs (depends on the MetricGraph package)
* Added cgeneric versions of the stationary models
* Replaced rgeneric models by their cgeneric counterparts
* Added a new parameterization (range and std. dev)
* Created a new method gg_df to help posterior plotting in ggplot2

# rSPDE 1.2.0
* Added an inlabru interface
* Added "rational.order" and "rational.type" functions
* Added the BRASIL rational approximation
* Improved covariance-based operator objects
* Improved log-likelihood computation
* Created 2d folded Matern under different boundary conditions
* Implemented different boundary conditions for 1d folded Matern


# rSPDE 1.1.1
* Adjusts on donttest examples for CRAN

# rSPDE 1.1.0
* Minor typos on vignettes and man pages were corrected
* Some examples were changed to improve their numerical stability

# rSPDE 1.0.0
* Implementation of the covariance-based rational approximation for stationary Matérn models
* R-INLA implementation of the rational SPDE approach
* Added an introduction to rSPDE vignette
* The previous vignette was updated an became an operator-based rational approximation vignette
* Added a vignette for the R-INLA implementation of the SPDE approach
* Added a vignette to present the rational approximation using the rSPDE package
* Backward compatibility was maintained

# rSPDE 0.6.3
* Change to inline citations in the Vignette to avoid problems on CRAN

# rSPDE 0.6.2

# rSPDE 0.6.1
* Add rgdal as suggested package

# rSPDE 0.5.0
* Remove dependency on INLA for Vignette on CRAN 
* Update citation 
