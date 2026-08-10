# Package index

## rSPDE package

- [`rSPDE`](https://davidbolin.github.io/rSPDE/reference/rSPDE-package.md)
  [`rSPDE-package`](https://davidbolin.github.io/rSPDE/reference/rSPDE-package.md)
  : Rational approximations of fractional SPDEs.

## rSPDE models

- [`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.md)
  : Matern rSPDE model object for INLA
- [`rspde.metric_graph()`](https://davidbolin.github.io/rSPDE/reference/rspde.metric_graph.md)
  : Matern rSPDE model object for metric graphs in INLA
- [`rspde.spacetime()`](https://davidbolin.github.io/rSPDE/reference/rspde.spacetime.md)
  : Space-Time Random Fields via SPDE Approximation
- [`rspde.anistropic2d()`](https://davidbolin.github.io/rSPDE/reference/rspde.anistropic2d.md)
  : Rational approximations of stationary anisotropic Gaussian Matern
  random fields
- [`rspde.matern1d()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern1d.md)
  : Matern rSPDE model object for INLA
- [`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.md)
  : Rational approximations of stationary Gaussian Matern random fields
- [`matern2d.operators()`](https://davidbolin.github.io/rSPDE/reference/matern2d.operators.md)
  : Rational approximations of stationary anisotropic Gaussian Matern
  random fields
- [`spde.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.operators.md)
  : Rational approximations of non-stationary Gaussian SPDE Matern
  random fields
- [`fractional.operators()`](https://davidbolin.github.io/rSPDE/reference/fractional.operators.md)
  : Rational approximations of fractional operators
- [`matern.rational()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.md)
  : Rational approximation of the Matern fields on intervals and metric
  graphs
- [`spacetime.operators()`](https://davidbolin.github.io/rSPDE/reference/spacetime.operators.md)
  : Space-time random fields
- [`hybrid.spde()`](https://davidbolin.github.io/rSPDE/reference/hybrid.spde.md)
  : Hybrid non-stationary Whittle-Matern SPDE model with non-zero mean
- [`rspde.hybrid.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.hybrid.matern.md)
  : Hybrid Whittle-Matern SPDE model for INLA / inlabru (alpha = 2)

## Linear mixed-effects models

- [`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.md)
  : rSPDE linear mixed effects models

- [`predict(`*`<rspde_lme>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.rspde_lme.md)
  : Prediction of a mixed effects regression model on a metric graph.

- [`summary(`*`<rspde_lme>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rspde_lme.md)
  :

  Summary Method for `rspde_lme` Objects.

- [`glance(`*`<rspde_lme>`*`)`](https://davidbolin.github.io/rSPDE/reference/glance.rspde_lme.md)
  :

  Glance at an `rspde_lme` object

- [`augment(`*`<rspde_lme>`*`)`](https://davidbolin.github.io/rSPDE/reference/augment.rspde_lme.md)
  :

  Augment data with information from a `rspde_lme` object

- [`posterior_crossvalidation()`](https://davidbolin.github.io/rSPDE/reference/posterior_crossvalidation.md)
  :

  Posterior cross-validation for `rspde_lme` models

## Intrinsic models

- [`intrinsic.matern.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.matern.operators.md)
  : Covariance-based approximations of intrinsic fields
- [`intrinsic.operators()`](https://davidbolin.github.io/rSPDE/reference/intrinsic.operators.md)
  : Covariance-based approximations of intrinsic fields
- [`variogram.intrinsic.spde()`](https://davidbolin.github.io/rSPDE/reference/variogram.intrinsic.spde.md)
  : Variogram of intrinsic SPDE model
- [`rspde.intrinsic.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.intrinsic.md)
  : Intrinsic Matern rSPDE model object for INLA
- [`rspde.intrinsic()`](https://davidbolin.github.io/rSPDE/reference/rspde.intrinsic.md)
  : Rational approximations of fractional intrinsic fields
- [`simulate(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.intrinsicCBrSPDEobj.md)
  : Simulation of a fractional intrinsic SPDE using the covariance-based
  rational SPDE approximation

## Log-likelihood

- [`rSPDE.matern.loglike()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.matern.loglike.md)
  : Object-based log-likelihood function for latent Gaussian fractional
  SPDE model using the rational approximations
- [`rSPDE.loglike()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.loglike.md)
  : Object-based log-likelihood function for latent Gaussian fractional
  SPDE model
- [`spde.matern.loglike()`](https://davidbolin.github.io/rSPDE/reference/spde.matern.loglike.md)
  : Parameter-based log-likelihood for a latent Gaussian Matern SPDE
  model using a rational SPDE approximation
- [`rSPDE.construct.matern.loglike()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.construct.matern.loglike.md)
  : Constructor of Matern loglikelihood functions.
- [`construct.spde.matern.loglike()`](https://davidbolin.github.io/rSPDE/reference/construct.spde.matern.loglike.md)
  : Constructor of Matern loglikelihood functions for non-stationary
  models.

## Computation of precision matrices

- [`rspde.matern.precision()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.precision.md)
  : Precision matrix of the covariance-based rational approximation of
  stationary Gaussian Matern random fields
- [`rspde.matern.precision.integer()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.precision.integer.md)
  : Precision matrix of stationary Gaussian Matern random fields with
  integer covariance exponent

## Methods for rSPDE, CBrSPDE and spacetime objects

- [`predict(`*`<rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.rSPDEobj.md)
  : Prediction of a fractional SPDE using a rational SPDE approximation
- [`simulate(`*`<rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.rSPDEobj.md)
  : Simulation of a fractional SPDE using a rational SPDE approximation
- [`summary(`*`<rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj.md)
  [`print(`*`<summary.rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj.md)
  [`print(`*`<rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj.md)
  : Summarise rSPDE objects
- [`update(`*`<rSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.rSPDEobj.md)
  : Update parameters of rSPDEobj objects
- [`predict(`*`<CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.CBrSPDEobj.md)
  : Prediction of a fractional SPDE using the covariance-based rational
  SPDE approximation
- [`simulate(`*`<CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.CBrSPDEobj.md)
  : Simulation of a fractional SPDE using the covariance-based rational
  SPDE approximation
- [`summary(`*`<CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj.md)
  [`print(`*`<summary.CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj.md)
  [`print(`*`<CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj.md)
  : Summarise CBrSPDE objects
- [`update(`*`<CBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.CBrSPDEobj.md)
  : Update parameters of CBrSPDEobj objects
- [`precision()`](https://davidbolin.github.io/rSPDE/reference/precision.CBrSPDEobj.md)
  : Get the precision matrix of CBrSPDEobj objects
- [`predict(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.CBrSPDEobj2d.md)
  : Prediction of an anisotropic Whittle-Matern field
- [`simulate(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.CBrSPDEobj2d.md)
  : Simulation of a fractional SPDE using the covariance-based rational
  SPDE approximation
- [`summary(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj2d.md)
  [`print(`*`<summary.CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj2d.md)
  [`print(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.CBrSPDEobj2d.md)
  : Summarise CBrSPDEobj2d objects
- [`update(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.CBrSPDEobj2d.md)
  : Update parameters of CBrSPDEobj2d objects
- [`precision(`*`<CBrSPDEobj2d>`*`)`](https://davidbolin.github.io/rSPDE/reference/precision.CBrSPDEobj2d.md)
  : Get the precision matrix of CBrSPDEobj2d objects
- [`simulate(`*`<rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.rSPDEobj1d.md)
  : Simulation of a Matern field using a rational SPDE approximation
- [`summary(`*`<rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj1d.md)
  [`print(`*`<summary.rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj1d.md)
  [`print(`*`<rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rSPDEobj1d.md)
  : Summarise rSPDE objects without FEM
- [`update(`*`<rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.rSPDEobj1d.md)
  : Update parameters of rSPDEobj1d objects
- [`precision(`*`<rSPDEobj1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/precision.rSPDEobj1d.md)
  : Get the precision matrix of rSPDEobj1d objects
- [`precision(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/precision.spacetimeobj.md)
  : Get the precision matrix of spacetimeobj objects
- [`predict(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.spacetimeobj.md)
  : Prediction of a space-time SPDE
- [`simulate(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.spacetimeobj.md)
  : Simulation of space-time models
- [`update(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.spacetimeobj.md)
  : Update parameters of spacetimeobj objects
- [`summary(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.spacetimeobj.md)
  [`print(`*`<summary.spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.spacetimeobj.md)
  [`print(`*`<spacetimeobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.spacetimeobj.md)
  : Summarise spacetime objects
- [`precision(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/precision.intrinsicCBrSPDEobj.md)
  : Get the precision matrix of intrinsicCBrSPDEobj objects
- [`predict(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.intrinsicCBrSPDEobj.md)
  : Prediction of an intrinsic Whittle-Matern model
- [`summary(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.intrinsicCBrSPDEobj.md)
  [`print(`*`<summary.intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.intrinsicCBrSPDEobj.md)
  [`print(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.intrinsicCBrSPDEobj.md)
  : Summary method for class "intrinsicCBrSPDEobj"
- [`update(`*`<intrinsicCBrSPDEobj>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.intrinsicCBrSPDEobj.md)
  : Update parameters of intrinsicCBrSPDEobj objects
- [`predict(`*`<hybrid_spde>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.hybrid_spde.md)
  : Kriging prediction for a hybrid Whittle-Matern SPDE model
- [`simulate(`*`<hybrid_spde>`*`)`](https://davidbolin.github.io/rSPDE/reference/simulate.hybrid_spde.md)
  : Simulation of a hybrid Whittle-Matern SPDE model
- [`update(`*`<hybrid_spde>`*`)`](https://davidbolin.github.io/rSPDE/reference/update.hybrid_spde.md)
  : Update parameters of hybrid_spde objects

## Functions and methods for R-INLA rSPDE objects

- [`rspde.make.A()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.A.md)
  : Observation/prediction matrices for rSPDE models.

- [`spde.make.A()`](https://davidbolin.github.io/rSPDE/reference/spde.make.A.md)
  : Observation/prediction matrices for rSPDE models with integer
  smoothness.

- [`rspde.make.index()`](https://davidbolin.github.io/rSPDE/reference/rspde.make.index.md)
  : rSPDE model index vector generation

- [`graph_data_rspde()`](https://davidbolin.github.io/rSPDE/reference/graph_data_rspde.md)
  : Data extraction from metric graphs for 'rSPDE' models

- [`rspde.mesh.project()`](https://davidbolin.github.io/rSPDE/reference/rspde.mesh.project.md)
  [`rspde.mesh.projector()`](https://davidbolin.github.io/rSPDE/reference/rspde.mesh.project.md)
  :

  Calculate a lattice projection to/from an `inla.mesh` for rSPDE
  objects

- [`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.md)
  : rSPDE result extraction from INLA estimation results

- [`summary(`*`<rspde_result>`*`)`](https://davidbolin.github.io/rSPDE/reference/summary.rspde_result.md)
  :

  Summary for posteriors of field parameters for an `inla_rspde` model
  from a `rspde_result` object

- [`precision(`*`<inla_rspde>`*`)`](https://davidbolin.github.io/rSPDE/reference/precision.inla_rspde.md)
  :

  Get the precision matrix of `inla_rspde` objects

- [`gg_df()`](https://davidbolin.github.io/rSPDE/reference/gg_df.md) :
  Data frame for result objects from R-INLA fitted models to be used in
  ggplot2

- [`gg_df(`*`<rspde_result>`*`)`](https://davidbolin.github.io/rSPDE/reference/gg_df.rspde_result.md)
  : Data frame for rspde_result objects to be used in ggplot2

## Functions and methods for rSPDE interface for inlabru

- [`cross_validation()`](https://davidbolin.github.io/rSPDE/reference/cross_validation.md)
  : Perform cross-validation on a list of fitted models.
- [`predict(`*`<inla_rspde_matern1d>`*`)`](https://davidbolin.github.io/rSPDE/reference/predict.inla_rspde_matern1d.md)
  : Predict method for 'inlabru' stationary Matern 1d models
- [`ibm_n.bru_mapper_inla_rspde_fintrinsic()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`ibm_values.bru_mapper_inla_rspde_fintrinsic()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`ibm_jacobian.bru_mapper_inla_rspde_fintrinsic()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`bru_get_mapper.inla_rspde()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`ibm_n.bru_mapper_inla_rspde()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`ibm_values.bru_mapper_inla_rspde()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  [`ibm_jacobian.bru_mapper_inla_rspde()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde.md)
  : rSPDE inlabru mapper
- [`bru_get_mapper.inla_rspde_spacetime()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_spacetime.md)
  : rSPDE space time inlabru mapper
- [`bru_get_mapper.inla_rspde_anisotropic2d()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_anisotropic2d.md)
  : rSPDE anisotropic inlabru mapper
- [`bru_get_mapper.inla_rspde_matern1d()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_matern1d.md)
  [`ibm_n.bru_mapper_inla_rspde_matern1d()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_matern1d.md)
  [`ibm_values.bru_mapper_inla_rspde_matern1d()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_matern1d.md)
  [`ibm_jacobian.bru_mapper_inla_rspde_matern1d()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_matern1d.md)
  : rSPDE stationary inlabru mapper
- [`bru_get_mapper.inla_rspde_fintrinsic()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.inla_rspde_fintrinsic.md)
  : rSPDE inlabru mapper
- [`bru_get_mapper.intrinsic_matern()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.intrinsic_matern.md)
  [`ibm_n.bru_mapper_intrinsic_matern()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.intrinsic_matern.md)
  [`ibm_values.bru_mapper_intrinsic_matern()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.intrinsic_matern.md)
  [`ibm_jacobian.bru_mapper_intrinsic_matern()`](https://davidbolin.github.io/rSPDE/reference/bru_get_mapper.intrinsic_matern.md)
  : rSPDE inlabru mapper

## Finite element-related functions

- [`rSPDE.A1d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.A1d.md)
  : Observation matrix for finite element discretization on R
- [`rSPDE.fem1d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.fem1d.md)
  : Finite element calculations for problems on R
- [`rSPDE.fem2d()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.fem2d.md)
  : Finite element calculations for problems in 2D
- [`rSPDE.Ast()`](https://davidbolin.github.io/rSPDE/reference/rSPDE.Ast.md)
  : Observation matrix for space-time models

## Auxiliary functions

- [`create_train_test_indices()`](https://davidbolin.github.io/rSPDE/reference/create_train_test_indices.md)
  : Create train and test splits for cross-validation
- [`get.initial.values.rSPDE()`](https://davidbolin.github.io/rSPDE/reference/get.initial.values.rSPDE.md)
  : Initial values for log-likelihood optimization in rSPDE models with
  a latent stationary Gaussian Matern model
- [`require.nowarnings()`](https://davidbolin.github.io/rSPDE/reference/require.nowarnings.md)
  : Warnings free loading of add-on packages
- [`matern.covariance()`](https://davidbolin.github.io/rSPDE/reference/matern.covariance.md)
  : The Matern covariance function
- [`matern.rational.cov()`](https://davidbolin.github.io/rSPDE/reference/matern.rational.cov.md)
  : Rational approximation of the Matern covariance
- [`folded.matern.covariance.1d()`](https://davidbolin.github.io/rSPDE/reference/folded.matern.covariance.1d.md)
  : The 1d folded Matern covariance function
- [`folded.matern.covariance.2d()`](https://davidbolin.github.io/rSPDE/reference/folded.matern.covariance.2d.md)
  : The 2d folded Matern covariance function
- [`rspde.matern.precision.opt()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.precision.opt.md)
  : Optimized precision matrix of the covariance-based rational
  approximation
- [`rspde.matern.precision.integer.opt()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.precision.integer.opt.md)
  : Optimized precision matrix of stationary Gaussian Matern random
  fields with integer covariance exponent
- [`` `rational.order<-`() ``](https://davidbolin.github.io/rSPDE/reference/rational.order-set.md)
  : Changing the order of the rational approximation
- [`` `rational.type<-`() ``](https://davidbolin.github.io/rSPDE/reference/rational.type-set.md)
  : Changing the type of the rational approximation
- [`rational.order()`](https://davidbolin.github.io/rSPDE/reference/rational.order.md)
  : Get the order of rational approximation.
- [`rational.type()`](https://davidbolin.github.io/rSPDE/reference/rational.type.md)
  : Get type of rational approximation.
- [`transform_parameters_anisotropic()`](https://davidbolin.github.io/rSPDE/reference/transform_parameters_anisotropic.md)
  : Transform Anisotropic SPDE Model Parameters to Original Scale
- [`transform_parameters_spacetime()`](https://davidbolin.github.io/rSPDE/reference/transform_parameters_spacetime.md)
  : Transform Spacetime SPDE Model Parameters to Original Scale
- [`cov_function_mesh()`](https://davidbolin.github.io/rSPDE/reference/cov_function_mesh.md)
  : Covariance between mesh nodes and locations
- [`covariance_mesh()`](https://davidbolin.github.io/rSPDE/reference/covariance_mesh.md)
  : Covariance between mesh nodes
- [`make_A()`](https://davidbolin.github.io/rSPDE/reference/make_A.md) :
  Projection matrix for model objects

## Operator operations

- [`Pr.mult()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Pr.solve()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Pl.mult()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Pl.solve()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Q.mult()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Q.solve()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Qsqrt.mult()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Qsqrt.solve()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Sigma.mult()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  [`Sigma.solve()`](https://davidbolin.github.io/rSPDE/reference/operator.operations.md)
  : Operations with the Pr and Pl operators
