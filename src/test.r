prepare_cgeneric_data <- function(data) {
  # Preprocess integers
  if (!is.null(data$ints)) {
    data$ints <- lapply(data$ints, function(int_vec) {
      if (length(int_vec) == 0 || is.null(int_vec)) {
        as.integer(NA) # Replace empty or NULL values with NA
      } else {
        as.integer(int_vec)
      }
    })
  } else {
    data$ints <- list() # Ensure ints is at least an empty list
  }

  # Preprocess doubles
  if (!is.null(data$doubles)) {
    data$doubles <- lapply(data$doubles, function(double_vec) {
      if (length(double_vec) == 0 || is.null(double_vec)) {
        as.numeric(NA) # Replace empty or NULL values with NA
      } else {
        as.numeric(double_vec)
      }
    })
  } else {
    data$doubles <- list() # Ensure doubles is at least an empty list
  }

  # Preprocess sparse matrices
  if (!is.null(data$smatrices)) {
    data$smatrices <- lapply(data$smatrices, function(sm) {
      nrow <- sm[1]
      ncol <- sm[2]
      n <- sm[3]
      if (n == 0) {
        return(NULL) # Handle completely empty sparse matrices
      }
      i <- sm[4:(4 + n - 1)]
      j <- sm[(4 + n):(4 + 2 * n - 1)]
      x <- sm[(4 + 2 * n):(4 + 3 * n - 1)]
      list(
        nrow = as.integer(nrow),
        ncol = as.integer(ncol),
        n = as.integer(n),
        i = as.integer(i),
        j = as.integer(j),
        x = as.numeric(x)
      )
    })
  } else {
    data$smatrices <- list() # Ensure smatrices is at least an empty list
  }

  # Preprocess dense matrices
  if (!is.null(data$matrices)) {
    data$matrices <- lapply(data$matrices, function(mat) {
      if (length(mat) < 2) {
        return(NULL) # Handle completely empty matrices
      }
      nrow <- mat[1]
      ncol <- mat[2]
      x <- mat[-(1:2)] # Exclude the first two elements (nrow, ncol)
      list(
        nrow = as.integer(nrow),
        ncol = as.integer(ncol),
        x = as.numeric(x)
      )
    })
  } else {
    data$matrices <- list() # Ensure matrices is at least an empty list
  }

  # Ensure characters is a character vector
  if (!is.character(data$characters)) {
    data$characters <- unlist(data$characters)
  }

  return(data)
}

call_dynamic_cgeneric_model <- function(model, cmd, theta) {
  # Validate model structure
  if (!is.list(model) || is.null(model$f) || is.null(model$f$cgeneric)) {
    stop("The model must be a list containing an $f$cgeneric structure.")
  }

  cgeneric <- model$f$cgeneric
  if (!is.list(cgeneric) || is.null(cgeneric$shlib) || is.null(cgeneric$model) || is.null(cgeneric$data)) {
    stop("The model$f$cgeneric structure must contain 'shlib', 'model', and 'data' elements.")
  }

  # Load the shared library explicitly
  so_path <- cgeneric$shlib
  if (!file.exists(so_path)) {
    stop(sprintf("Shared library '%s' not found.", so_path))
  }

  tryCatch({
    dyn.load(so_path)
  }, error = function(e) {
    stop(sprintf("Failed to load shared library '%s': %s", so_path, e$message))
  })

  # Preprocess data
  cgeneric$data <- prepare_cgeneric_data(cgeneric$data)

  # Convert cmd to integer
  cmd <- switch(
    cmd,
    "void" = 0,
    "Q" = 1,
    "graph" = 2,
    "mu" = 3,
    "initial" = 4,
    "log_norm_const" = 5,
    "log_prior" = 6,
    "quit" = 7,
    stop(sprintf("Invalid command '%s'", cmd))
  )

  # Validate theta
  if (!is.numeric(theta) || length(theta) == 0) {
    stop("The 'theta' argument must be a non-empty numeric vector.")
  }

  # Call the C function
  result <- tryCatch(
    {
      .Call(
        "call_dynamic_inla_cgeneric", # Call the C wrapper function
        as.integer(cmd),
        as.numeric(theta),
        cgeneric$data,
        cgeneric$model,
        cgeneric$shlib
      )
    },
    error = function(e) {
      stop(
        sprintf(
          "Error calling function '%s' in shared library '%s': %s",
          cgeneric$model,
          cgeneric$shlib,
          e$message
        )
      )
    }
  )

  # Optionally unload the shared library after use
  tryCatch({
    dyn.unload(so_path)
  }, error = function(e) {
    warning(sprintf("Failed to unload shared library '%s': %s", so_path, e$message))
  })

  return(result)
}

library(fmesher)
n_loc <- 20
loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
mesh_2d <- fm_mesh_2d(
    loc = loc_2d_mesh,
    cutoff = 0.01,
    max.edge = c(0.1, 0.5)
)
library(inlabru)
library(rSPDE)
model_aniso <- rspde.anistropic2d(mesh = mesh_2d)
model_1 <- rspde.matern(mesh = mesh_2d)


call_dynamic_cgeneric_model(model_1, cmd = "Q", theta = c(-1,-1,0,0,0))
