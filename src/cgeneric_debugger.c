#include <R.h>
#include <Rinternals.h>
#include <dlfcn.h>
#include <stdio.h>
#include "cgeneric.h"

SEXP call_dynamic_inla_cgeneric(SEXP r_cmd, SEXP r_theta, SEXP r_data, SEXP r_func_name, SEXP r_so_path) {
    // Validate inputs
    if (!isString(r_func_name) || LENGTH(r_func_name) != 1) {
        Rf_error("Function name must be a single string.");
    }
    if (!isString(r_so_path) || LENGTH(r_so_path) != 1) {
        Rf_error("Shared library path must be a single string.");
    }

    const char *func_name = CHAR(STRING_ELT(r_func_name, 0));
    const char *so_path = CHAR(STRING_ELT(r_so_path, 0));

    // Dynamically load the shared library
    void *handle = dlopen(so_path, RTLD_LAZY);
    if (!handle) {
        Rf_error("Failed to load shared library '%s': %s", so_path, dlerror());
    }

    // Lookup the function
    inla_cgeneric_func_tp *dynamic_function;
    *(void **)(&dynamic_function) = dlsym(handle, func_name);
    if (!dynamic_function) {
        dlclose(handle);
        Rf_error("Function '%s' not found in shared library '%s': %s", func_name, so_path, dlerror());
    }

    // Parse cmd
    if (!isInteger(r_cmd)) Rf_error("Expected integer for 'cmd'.");
    int cmd = INTEGER(r_cmd)[0];
    Rprintf("Received cmd = %d\n", cmd);
    if (cmd < INLA_CGENERIC_VOID || cmd > INLA_CGENERIC_QUIT) {
        dlclose(handle);
        Rf_error("Invalid cmd value: %d. Must be within enum inla_cgeneric_cmd_tp.", cmd);
    }
    inla_cgeneric_cmd_tp cmd_parsed = (inla_cgeneric_cmd_tp)cmd;

    // Print the parsed command name
    const char *cmd_name = INLA_CGENERIC_CMD_NAME(cmd_parsed);
    Rprintf("Parsed cmd = %d, Command Name = '%s'\n", cmd_parsed, cmd_name);

    // Parse theta
    if (!isReal(r_theta)) {
        Rf_error("Expected numeric vector for 'theta'.");
    }
    double *theta = REAL(r_theta);
    int theta_len = LENGTH(r_theta);
    Rprintf("Received theta of length %d\n", theta_len);
    
    // Print values of theta
    Rprintf("Theta values: ");
    for (int i = 0; i < theta_len; i++) {
        Rprintf("%f ", theta[i]);
    }
    Rprintf("\n");

    // Debug r_data structure
    Rprintf("Debugging r_data structure:\n");
    for (int i = 0; i < LENGTH(r_data); i++) {
        SEXP current_element = VECTOR_ELT(r_data, i);
        Rprintf("Element %d: TYPEOF = %d, LENGTH = %d\n", i, TYPEOF(current_element), LENGTH(current_element));
    }

    inla_cgeneric_data_tp data = {0};

    // Parse integers
    SEXP r_ints = VECTOR_ELT(r_data, 0);
    data.n_ints = LENGTH(r_ints);
    Rprintf("r_ints length: %d\n", data.n_ints);

    if (data.n_ints > 0) {
        data.ints = Calloc(data.n_ints, inla_cgeneric_vec_tp *);
        for (int i = 0; i < data.n_ints; i++) {
            SEXP current_int = VECTOR_ELT(r_ints, i);
            SEXP name = getAttrib(r_ints, R_NamesSymbol);
            const char *int_name = (name != R_NilValue) ? CHAR(STRING_ELT(name, i)) : "unnamed";
            if (TYPEOF(current_int) != INTSXP) {
                Rf_error("ints[%d] (%s) is not of type integer (TYPEOF = %d).", i, int_name, TYPEOF(current_int));
            }
            data.ints[i] = Calloc(1, inla_cgeneric_vec_tp);
            data.ints[i]->len = LENGTH(current_int);
            data.ints[i]->ints = INTEGER(current_int);
            data.ints[i]->name = NULL;
            Rprintf("Successfully parsed ints[%d] (%s):\n", i, int_name);
            for (int j = 0; j < data.ints[i]->len && j < 10; j++) {
                Rprintf("  %d\n", data.ints[i]->ints[j]);
            }
        }
    } else {
        data.ints = NULL; // No integers present
    }

    // Parse doubles
    SEXP r_doubles = VECTOR_ELT(r_data, 1);
    data.n_doubles = LENGTH(r_doubles);
    Rprintf("r_doubles length: %d\n", data.n_doubles);

    if (data.n_doubles > 0) {
        data.doubles = Calloc(data.n_doubles, inla_cgeneric_vec_tp *);
        for (int i = 0; i < data.n_doubles; i++) {
            SEXP current_double = VECTOR_ELT(r_doubles, i);
            SEXP name = getAttrib(r_doubles, R_NamesSymbol);
            const char *double_name = (name != R_NilValue) ? CHAR(STRING_ELT(name, i)) : "unnamed";
            if (TYPEOF(current_double) != REALSXP) {
                Rf_error("doubles[%d] (%s) is not of type numeric (TYPEOF = %d).", i, double_name, TYPEOF(current_double));
            }
            data.doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
            data.doubles[i]->len = LENGTH(current_double);
            data.doubles[i]->doubles = REAL(current_double);
            data.doubles[i]->name = NULL;
            Rprintf("Successfully parsed doubles[%d] (%s):\n", i, double_name);
            for (int j = 0; j < data.doubles[i]->len && j < 10; j++) {
                Rprintf("  %f\n", data.doubles[i]->doubles[j]);
            }
        }
    } else {
        data.doubles = NULL; // No doubles present
    }

    // Parse dense matrices
    SEXP r_mats = VECTOR_ELT(r_data, 3);
    data.n_mats = LENGTH(r_mats);
    Rprintf("r_mats length: %d\n", data.n_mats);

    if (data.n_mats > 0) {
        data.mats = Calloc(data.n_mats, inla_cgeneric_mat_tp *);
        for (int i = 0; i < data.n_mats; i++) {
            SEXP current_mat = VECTOR_ELT(r_mats, i);
            SEXP name = getAttrib(r_mats, R_NamesSymbol);
            const char *mat_name = (name != R_NilValue) ? CHAR(STRING_ELT(name, i)) : "unnamed";
            Rprintf("mats[%d] (%s): TYPEOF = %d, LENGTH = %d\n", i, mat_name, TYPEOF(current_mat), LENGTH(current_mat));

            if (current_mat == R_NilValue) {
                Rprintf("mats[%d] (%s) is NULL\n", i, mat_name);
                continue;
            }

            // Dense matrices are VECSXP and must be handled as structured objects
            if (TYPEOF(current_mat) != VECSXP) {
                Rf_error("mats[%d] (%s) is not a structured object (TYPEOF = %d).", i, mat_name, TYPEOF(current_mat));
            }

            // Extract nrow, ncol, and x
            SEXP r_nrow = VECTOR_ELT(current_mat, 0);
            SEXP r_ncol = VECTOR_ELT(current_mat, 1);
            SEXP r_x = VECTOR_ELT(current_mat, 2);

            if (TYPEOF(r_nrow) != INTSXP || LENGTH(r_nrow) != 1) {
                Rf_error("mats[%d] (%s) has invalid 'nrow' (TYPEOF = %d, LENGTH = %d).", i, mat_name, TYPEOF(r_nrow), LENGTH(r_nrow));
            }
            if (TYPEOF(r_ncol) != INTSXP || LENGTH(r_ncol) != 1) {
                Rf_error("mats[%d] (%s) has invalid 'ncol' (TYPEOF = %d, LENGTH = %d).", i, mat_name, TYPEOF(r_ncol), LENGTH(r_ncol));
            }
            if (TYPEOF(r_x) != REALSXP) {
                Rf_error("mats[%d] (%s) has invalid 'x' (TYPEOF = %d).", i, mat_name, TYPEOF(r_x));
            }

            // Allocate and parse the matrix
            data.mats[i] = Calloc(1, inla_cgeneric_mat_tp);
            data.mats[i]->nrow = INTEGER(r_nrow)[0];
            data.mats[i]->ncol = INTEGER(r_ncol)[0];
            data.mats[i]->x = REAL(r_x);

            Rprintf("Successfully parsed mats[%d] (%s): nrow=%d, ncol=%d\n", i, mat_name, data.mats[i]->nrow, data.mats[i]->ncol);
            for (int j = 0; j < data.mats[i]->nrow * data.mats[i]->ncol && j < 10; j++) {
                Rprintf("  x[%d]=%f\n", j, data.mats[i]->x[j]);
            }
        }
    } else {
        data.mats = NULL; // No dense matrices present
    }
    // Parse sparse matrices
    SEXP r_smats = VECTOR_ELT(r_data, 4);
    data.n_smats = LENGTH(r_smats);
    Rprintf("r_smats length: %d\n", data.n_smats);

    if (data.n_smats > 0) {
        data.smats = Calloc(data.n_smats, inla_cgeneric_smat_tp *);
        for (int i = 0; i < data.n_smats; i++) {
            SEXP current_smat = VECTOR_ELT(r_smats, i);
            SEXP name = getAttrib(r_smats, R_NamesSymbol);
            const char *smat_name = (name != R_NilValue) ? CHAR(STRING_ELT(name, i)) : "unnamed";
            Rprintf("smats[%d] (%s): TYPEOF = %d, LENGTH = %d\n", i, smat_name, TYPEOF(current_smat), LENGTH(current_smat));

            if (current_smat == R_NilValue) {
                Rprintf("smats[%d] (%s) is NULL\n", i, smat_name);
                continue;
            }

            // Sparse matrices are structured as lists (VECSXP)
            if (TYPEOF(current_smat) != VECSXP || LENGTH(current_smat) != 6) {
                Rf_error("smats[%d] (%s) does not have the required 6 components (TYPEOF = %d, LENGTH = %d).",
                         i, smat_name, TYPEOF(current_smat), LENGTH(current_smat));
            }

            // Extract components: nrow, ncol, n, i, j, and x
            SEXP r_nrow = VECTOR_ELT(current_smat, 0);
            SEXP r_ncol = VECTOR_ELT(current_smat, 1);
            SEXP r_n = VECTOR_ELT(current_smat, 2);
            SEXP r_i = VECTOR_ELT(current_smat, 3);
            SEXP r_j = VECTOR_ELT(current_smat, 4);
            SEXP r_x = VECTOR_ELT(current_smat, 5);

            // Debugging component types and lengths
            Rprintf("smats[%d] (%s) components:\n", i, smat_name);
            Rprintf("  nrow: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_nrow), LENGTH(r_nrow));
            Rprintf("  ncol: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_ncol), LENGTH(r_ncol));
            Rprintf("  n: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_n), LENGTH(r_n));
            Rprintf("  i: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_i), LENGTH(r_i));
            Rprintf("  j: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_j), LENGTH(r_j));
            Rprintf("  x: TYPEOF = %d, LENGTH = %d\n", TYPEOF(r_x), LENGTH(r_x));

            // Validate types
            if (TYPEOF(r_nrow) != INTSXP || LENGTH(r_nrow) != 1) {
                Rf_error("smats[%d] (%s) has invalid 'nrow' (TYPEOF = %d, LENGTH = %d).", i, smat_name, TYPEOF(r_nrow), LENGTH(r_nrow));
            }
            if (TYPEOF(r_ncol) != INTSXP || LENGTH(r_ncol) != 1) {
                Rf_error("smats[%d] (%s) has invalid 'ncol' (TYPEOF = %d, LENGTH = %d).", i, smat_name, TYPEOF(r_ncol), LENGTH(r_ncol));
            }
            if (TYPEOF(r_n) != INTSXP || LENGTH(r_n) != 1) {
                Rf_error("smats[%d] (%s) has invalid 'n' (TYPEOF = %d, LENGTH = %d).", i, smat_name, TYPEOF(r_n), LENGTH(r_n));
            }
            if (TYPEOF(r_i) != INTSXP) {
                Rf_error("smats[%d] (%s) has invalid 'i' (TYPEOF = %d).", i, smat_name, TYPEOF(r_i));
            }
            if (TYPEOF(r_j) != INTSXP) {
                Rf_error("smats[%d] (%s) has invalid 'j' (TYPEOF = %d).", i, smat_name, TYPEOF(r_j));
            }
            if (TYPEOF(r_x) != REALSXP) {
                Rf_error("smats[%d] (%s) has invalid 'x' (TYPEOF = %d).", i, smat_name, TYPEOF(r_x));
            }

            // Allocate and parse the sparse matrix
            data.smats[i] = Calloc(1, inla_cgeneric_smat_tp);
            data.smats[i]->nrow = INTEGER(r_nrow)[0];
            data.smats[i]->ncol = INTEGER(r_ncol)[0];
            data.smats[i]->n = INTEGER(r_n)[0];
            data.smats[i]->i = INTEGER(r_i);
            data.smats[i]->j = INTEGER(r_j);
            data.smats[i]->x = REAL(r_x);

            // Debug parsed values
            Rprintf("Successfully parsed smats[%d] (%s): nrow=%d, ncol=%d, n=%d\n",
                    i, smat_name, data.smats[i]->nrow, data.smats[i]->ncol, data.smats[i]->n);
            for (int j = 0; j < data.smats[i]->n && j < 10; j++) {
                Rprintf("  i[%d]=%d, j[%d]=%d, x[%d]=%f\n",
                        j, data.smats[i]->i[j], j, data.smats[i]->j[j], j, data.smats[i]->x[j]);
            }
        }
    } else {
        data.smats = NULL; // No sparse matrices present
    }

    Rprintf("All data successfully parsed. Calling the dynamic function now...\n");

    // Dynamically call the loaded function
    double *result = dynamic_function(cmd_parsed, theta, &data);

    // Convert the result to an R vector
    SEXP r_result = PROTECT(allocVector(REALSXP, theta_len));
    for (int i = 0; i < theta_len; i++) {
        REAL(r_result)[i] = result[i];
    }
    UNPROTECT(1);

    // Free allocated memory
    if (data.ints) {
        for (int i = 0; i < data.n_ints; i++) {
            Free(data.ints[i]);
        }
        Free(data.ints);
    }
    if (data.doubles) {
        for (int i = 0; i < data.n_doubles; i++) {
            Free(data.doubles[i]);
        }
        Free(data.doubles);
    }
    if (data.mats) {
        for (int i = 0; i < data.n_mats; i++) {
            Free(data.mats[i]);
        }
        Free(data.mats);
    }
    if (data.smats) {
        for (int i = 0; i < data.n_smats; i++) {
            Free(data.smats[i]);
        }
        Free(data.smats);
    }

    // Close the shared library
    dlclose(handle);

    return r_result;
}