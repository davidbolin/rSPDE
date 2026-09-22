# rSPDE (development version)

* Added `type_rational_approximation = "wl2"`, a new way of obtaining the
  rational coefficients. The classes are parameterised so that every fit is a 
  valid model and have no constant term, so the covariance-based models have `m` 
  instead of `m + 1` latent blocks. Available in `matern.operators()` and
  `CBrSPDE.matern.operators()` for both `type = "covariance"` and
  `type = "operator"`, in `matern.rational()` and `matern.rational.cov()` for
  the exact one-dimensional models, and directly through
  `rational.coefficients.wl2()`. It requires `floor(nu + d/2)` to be 0 or 1 for
  the covariance type and `nu + d/2 < 2` for the operator type. The existing
  `"brasil"`, `"chebfun"` and `"chebfunLB"` types are unchanged and remain the
  default.
* The mesh-free weighted-L2 coefficients are stored in the package, as the
  tabulated ones are, for `d` = 1 to 3, `m` = 1 to 6 and `floor(alpha)` = 0 and
  1. They depend on neither the mesh nor `kappa` and are what a covariance-based 
  model uses by default, so the common case fits nothing at set-up. 
  `data-raw/wl2_tables.R` regenerates them, and the tests check that a fresh fit 
  still reproduces what is stored.
* Added `rspde.xmin()`, which computes the lower end of the spectral interval
  from the mesh and a lower bound for `kappa`, and the arguments `x_min` and
  `kappa_ref` to `matern.operators()`. Fitting on that shorter interval is
  appreciably more accurate than the mesh-free fit, and is done when the model
  is created. `update_rational_coefficients()` recomputes the coefficients once
  `kappa` has been estimated.
* Added `rspde.wl2.table()`, which builds a weighted-L2 table for a given
  spectral interval, and the `wl2_table` argument of `matern.operators()`,
  which takes one and overrides `x_min` and `kappa_ref`. A table built for
  another dimension, order or range of `alpha` is refused rather than used.
* Added `rspde.cache()`, which keeps generated coefficient tables between
  sessions, under `tools::R_user_dir("rSPDE", "cache")` or a directory of your
  choosing. It is off by default, since a package should not write outside the
  session temporary directory unless asked; the environment variable
  `RSPDE_CACHE_DIR` sets it for non-interactive use. Lookup is automatic, and a
  table fitted on a slightly wider spectral interval is reused, since it still
  covers the whole spectrum.
* Added `variance_correction = "nodal"` to `matern.operators()`, which adds
  `max(sigma^2 - diag(Sigma), 0)` to the diagonal of the covariance of a
  `"wl2"` covariance-based model. What this corrects is mostly the finite
  element discretisation rather than the rational approximation. It is off by 
  default, as it depends on the parameters it cannot be tabulated, so it is meant 
  for a model whose parameters are already estimated.
* The package now has compiled code in every install, CRAN included:
  `src/wl2_fit.cpp` holds the inner loop of the weighted-L2 fit. The INLA
  `cgeneric` sources remain optional and are still built only with
  `RSPDE_COMPILE=1` or `--configure-args='--enable-compiled'`. The equivalent R
  implementation is kept as the reference and is used when
  `options(rSPDE.wl2.use.cpp = FALSE)`.
* The `RSpectra` dependency is gone. `rspde.xmin(eigenvalue = "exact")` and the
  scaling of the intrinsic models now use Lanczos iterations in the package.
  Both are also more robust, and the intrinsic one is faster than the old method.
* `intrinsic.operators()` now honours its `opts` argument, which was built and
  then replaced by a hardcoded list. Its entries are `tol` and `maxitr`, as
  before.
* Fixed `matern.rational.cov()`, which evaluated the covariance at the lags
  `h[1] - h`, rather than at `h`. A matrix of lags is now also accepted, and 
  returns a matrix.
* `get.roots()` now uses spline interpolation by default. Linear interpolation
  lost accuracy off the 200-node beta grid of the tables (symbol error 1.6e-4
  instead of 1.8e-6 at beta = 0.875, m = 4).
* Fixed the inlabru mapper for `rspde.spacetime()` models whose spatial mesh is
  a `metric_graph`. `bru_get_mapper()` used `bm_fmesher()` for the graph, which
  has no `fm_dof()` method, so `bru()` failed with
  "invalid subscript type 'list'".
* The optional cgeneric `Makefile` now only uses Homebrew `gcc-14` on macOS and
  the compilers R was configured with elsewhere, so compiled installs
  (`RSPDE_COMPILE=1`) work on Linux.
* `posterior_crossvalidation()` is now an S3 generic, with methods for
  `rspde_lme` fits and for lists of fitted models. MetricGraph provides the
  `graph_lme` method, so the two packages no longer mask each other's
  function, and a list can mix `rspde_lme` and `graph_lme` fits.

# rSPDE 2.6.0

* Added `posterior_crossvalidation()` for objects fitted with `rspde_lme()`.
  The function mirrors the interface of `MetricGraph::posterior_crossvalidation`.
* Added `hybrid.spde()`, a new hybrid Whittle-Matern SPDE model with a
  non-zero deterministic mean. 
* Added `rspde.hybrid.matern()`, a  INLA cgeneric model for the
  hybrid Whittle-Matern SPDE with alpha = 2.
* Added a `kappa_mu` option that lets the operator applied to the mean
  use a different range parameter from the one in the covariance.
  Default (`kappa_mu = NULL` in `hybrid.spde()`, `separate_kappa_mu =
  FALSE` in `rspde.hybrid.matern()`) keeps them linked. When enabled,
  `kappa_mu` is estimated jointly in `rspde_lme` and INLA, or can be
  held fixed via `model_options$fix_kappa_mu`. 
* The INLA cgeneric models now use the rSPDE models built into INLA when they
  are available, falling back to the local rSPDE shared library otherwise.
  `shared_lib` also accepts a path to a shared library file.
* Added `rspde_safe_inla()` and `local_rspde_safe_inla()`, which check that a
  usable INLA installation is available, for use in examples and tests.
* `rspde.metric_graph()` now passes `shared_lib` on to `rspde.matern()`, and
  its default is now `"detect"`, matching the other INLA models.
* Updated the inlabru interface to the inlabru 2.14 API. rSPDE now requires
  `fmesher (>= 0.7.0)` and suggests `inlabru (>= 2.14.0)`.
* `predict.rspde_lme()` can reuse precomputed parameter-dependent quantities,
  which speeds up repeated predictions such as in cross-validation.
* Fixed `predict.rspde_lme()` for models with replicates: the same location in
  different replicates no longer triggers a duplicated-locations warning.
* Fixed `update()` for non-stationary models, where new `theta` values were
  ignored, so predictions from non-stationary `rspde_lme()` fits used stale
  parameters.
* `spde.matern.operators()` no longer converts a model to a stationary one when
  `B.tau` or `B.kappa` vary in space.
* Added a vignette comparing rSPDE with the exact Matern covariance in terms of
  timing and memory.

# rSPDE 2.5.2

* Added intrinsic Matern mapper support and related documentation.
* Added `covariance_mesh`, `cov_function_mesh`, and `make_A` documentation.
* Improved intrinsic and fractional operator implementations and stability.
* Updated INLA/inlabru interfaces and examples.
* Expanded unit tests and vignette updates (spacetime/anisotropic/intrinsic).

# rSPDE 2.5.1

* Added `model_options` argument to `rspde_lme()` function, which allows users to set starting values for different parameters and also to fix parameters during estimation.
* Added `previous_fit` argument to `rspde_lme()`, which allows users to provide a previously fitted model as input to obtain starting values for a new fit.

# rSPDE 2.5.0

* Improved the `cross_validation` function to allow for multiple likelihoods. 
* General adjusts on `rspde.intrinsic` for stability. 
* inlabru implementation for `rspde.intrinsic`. 
* Improved warning messages when calling inla-related functions. 
* Added `wCRPS` and `swCRPS` scores on `cross_validation`.

# rSPDE 2.4.0

* Created the `group_predict` function, to obtain predictions on a testing set based on observations on a training set.
* Added support for `stochvol`, `stochvol.nig`, `stochvolln` and `binomial` likelihoods in `cross_validation` function.
* Changing the default `nu.upper.bound` to 2 in dimension 1, and keeping the default `nu.upper.bound` to 4 in dimension 2 in `rspde.matern()` function.
* Created `matern.rational()` operators for creating stationary matern operators.
* Created `spacetime.operators()` for creating space-time models.
* Created `matern2d.operators()` for anisotropic operators.
* Implemented space-time operators in cgeneric to be used in `INLA` and `inlabru`.
* Implemented anisotropic operators in cgeneric to be used in `INLA` and `inlabru`.
* Implemented stationary operators in cgeneric to be used in `INLA` and `inlabru`.
* Added vignette on space-time models.
* Added vignette on stationary models.
* Added vignette on anisotropic models.

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
