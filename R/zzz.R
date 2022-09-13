# from: https://github.com/tidyverse/hms/blob/master/R/zzz.R
# Thu Apr 19 10:53:24 CEST 2018
register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1)
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }

  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

register_all_s3_methods = function() {
  inlabru_installed <- "inlabru" %in% rownames(installed.packages())
  if(inlabru_installed){
    register_s3_method("inlabru", "bru_mapper", "inla_rspde")
    register_s3_method("inlabru", "ibm_n", "bru_mapper_inla_rspde") 
    register_s3_method("inlabru", "ibm_values", "bru_mapper_inla_rspde")
    register_s3_method("inlabru", "ibm_amatrix", "bru_mapper_inla_rspde")

    inlabru_version <- packageVersion("inlabru")
    if(inlabru_version >= "2.5.3.9002"){
    register_s3_method("inlabru", "bru_get_mapper", "inla_rspde")
    }
  }
}

.onLoad = function(libname, pkgname) {
  register_all_s3_methods()
}