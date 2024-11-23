#include <R.h>
#include <Rinternals.h>
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
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

    // Parse theta
    if (!isReal(r_theta)) Rf_error("Expected numeric vector for 'theta'.");
    double *theta = REAL(r_theta);
    int theta_len = LENGTH(r_theta);
    Rprintf("Received theta of length %d\n", theta_len);

    inla_cgeneric_data_tp data = {0};

    // Parse integers
    SEXP r_ints = VECTOR_ELT(r_data, 0);
    data.n_ints = LENGTH(r_ints);
    Rprintf("r_ints length: %d\n", data.n_ints);

    if (data.n_ints > 0) {
        data.ints = Calloc(data.n_ints, inla_cgeneric_vec_tp *);
        SEXP int_names = getAttrib(r_ints, R_NamesSymbol);
        for (int i = 0; i < data.n_ints; i++) {
            SEXP current_int = VECTOR_ELT(r_ints, i);
            if (TYPEOF(current_int) != INTSXP) {
                Rf_error("ints[%d] is not of type integer.", i);
            }
            data.ints[i] = Calloc(1, inla_cgeneric_vec_tp);
            data.ints[i]->len = LENGTH(current_int);
            data.ints[i]->ints = INTEGER(current_int);
            if (int_names != R_NilValue && LENGTH(int_names) > i) {
                data.ints[i]->name = strdup(CHAR(STRING_ELT(int_names, i)));
            } else {
                data.ints[i]->name = strdup("unnamed");
            }
            Rprintf("ints[%d] (%s): ", i, data.ints[i]->name);
            for (int j = 0; j < data.ints[i]->len; j++) {
                Rprintf("%d ", data.ints[i]->ints[j]);
            }
            Rprintf("\n");
        }
    }

    // Parse doubles
    SEXP r_doubles = VECTOR_ELT(r_data, 1);
    data.n_doubles = LENGTH(r_doubles);
    Rprintf("r_doubles length: %d\n", data.n_doubles);

    if (data.n_doubles > 0) {
        data.doubles = Calloc(data.n_doubles, inla_cgeneric_vec_tp *);
        SEXP double_names = getAttrib(r_doubles, R_NamesSymbol);
        for (int i = 0; i < data.n_doubles; i++) {
            SEXP current_double = VECTOR_ELT(r_doubles, i);
            if (TYPEOF(current_double) != REALSXP) {
                Rf_error("doubles[%d] is not of type numeric.", i);
            }
            data.doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
            data.doubles[i]->len = LENGTH(current_double);
            data.doubles[i]->doubles = REAL(current_double);
            if (double_names != R_NilValue && LENGTH(double_names) > i) {
                data.doubles[i]->name = strdup(CHAR(STRING_ELT(double_names, i)));
            } else {
                data.doubles[i]->name = strdup("unnamed");
            }
            Rprintf("doubles[%d] (%s): ", i, data.doubles[i]->name);
            for (int j = 0; j < data.doubles[i]->len; j++) {
                Rprintf("%f ", data.doubles[i]->doubles[j]);
            }
            Rprintf("\n");
        }
    }

    // Parse character arrays
    SEXP r_chars = VECTOR_ELT(r_data, 2); // Assuming character arrays are in the third slot
    data.n_chars = LENGTH(r_chars);
    Rprintf("r_chars length: %d\n", data.n_chars);

    if (data.n_chars > 0) {
        data.chars = Calloc(data.n_chars, inla_cgeneric_vec_tp *);
        SEXP char_names = getAttrib(r_chars, R_NamesSymbol);

        for (int i = 0; i < data.n_chars; i++) {
            SEXP current_char = VECTOR_ELT(r_chars, i);
            if (TYPEOF(current_char) != STRSXP || LENGTH(current_char) != 1) {
                Rf_error("chars[%d] must be a single string.", i);
            }

            data.chars[i] = Calloc(1, inla_cgeneric_vec_tp);
            data.chars[i]->len = strlen(CHAR(STRING_ELT(current_char, 0))); // Length of the string
            data.chars[i]->chars = strdup(CHAR(STRING_ELT(current_char, 0))); // Store the string

            // Store the name of the vector
            if (char_names != R_NilValue && LENGTH(char_names) > i) {
                data.chars[i]->name = strdup(CHAR(STRING_ELT(char_names, i)));
            } else {
                data.chars[i]->name = strdup("unnamed");
            }

            // Print debug info
            Rprintf("chars[%d] (%s): %s\n", i, data.chars[i]->name, data.chars[i]->chars);
        }
    } else {
        data.chars = NULL; // No character arrays present
    }

        // Parse dense matrices
    SEXP r_mats = VECTOR_ELT(r_data, 3);
    data.n_mats = LENGTH(r_mats);
    Rprintf("r_mats length: %d\n", data.n_mats);

    if (data.n_mats > 0) {
        data.mats = Calloc(data.n_mats, inla_cgeneric_mat_tp *);
        SEXP mat_names = getAttrib(r_mats, R_NamesSymbol);
        for (int i = 0; i < data.n_mats; i++) {
            SEXP current_mat = VECTOR_ELT(r_mats, i);
            if (TYPEOF(current_mat) != VECSXP || LENGTH(current_mat) != 3) {
                Rf_error("mats[%d] is not a valid dense matrix structure.", i);
            }
            SEXP r_nrow = VECTOR_ELT(current_mat, 0);
            SEXP r_ncol = VECTOR_ELT(current_mat, 1);
            SEXP r_x = VECTOR_ELT(current_mat, 2);

            if (TYPEOF(r_nrow) != INTSXP || LENGTH(r_nrow) != 1) {
                Rf_error("mats[%d] has invalid 'nrow'.", i);
            }
            if (TYPEOF(r_ncol) != INTSXP || LENGTH(r_ncol) != 1) {
                Rf_error("mats[%d] has invalid 'ncol'.", i);
            }
            if (TYPEOF(r_x) != REALSXP) {
                Rf_error("mats[%d] has invalid 'x'.", i);
            }

            data.mats[i] = Calloc(1, inla_cgeneric_mat_tp);
            data.mats[i]->nrow = INTEGER(r_nrow)[0];
            data.mats[i]->ncol = INTEGER(r_ncol)[0];
            data.mats[i]->x = REAL(r_x);

            if (mat_names != R_NilValue && LENGTH(mat_names) > i) {
                data.mats[i]->name = strdup(CHAR(STRING_ELT(mat_names, i)));
            } else {
                data.mats[i]->name = strdup("unnamed");
            }

            Rprintf("mats[%d] (%s): nrow=%d, ncol=%d\n", i, data.mats[i]->name, data.mats[i]->nrow, data.mats[i]->ncol);
            for (int j = 0; j < data.mats[i]->nrow * data.mats[i]->ncol && j < 10; j++) {
                Rprintf("  x[%d] = %f\n", j, data.mats[i]->x[j]);
            }
            if (data.mats[i]->nrow * data.mats[i]->ncol > 10) {
                Rprintf("  ... (only showing the first 10 elements)\n");
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
        SEXP smat_names = getAttrib(r_smats, R_NamesSymbol);
        for (int i = 0; i < data.n_smats; i++) {
            SEXP current_smat = VECTOR_ELT(r_smats, i);
            if (TYPEOF(current_smat) != VECSXP || LENGTH(current_smat) != 6) {
                Rf_error("smats[%d] is not a valid sparse matrix structure.", i);
            }
            SEXP r_nrow = VECTOR_ELT(current_smat, 0);
            SEXP r_ncol = VECTOR_ELT(current_smat, 1);
            SEXP r_n = VECTOR_ELT(current_smat, 2);
            SEXP r_i = VECTOR_ELT(current_smat, 3);
            SEXP r_j = VECTOR_ELT(current_smat, 4);
            SEXP r_x = VECTOR_ELT(current_smat, 5);

            if (TYPEOF(r_nrow) != INTSXP || LENGTH(r_nrow) != 1) {
                Rf_error("smats[%d] has invalid 'nrow'.", i);
            }
            if (TYPEOF(r_ncol) != INTSXP || LENGTH(r_ncol) != 1) {
                Rf_error("smats[%d] has invalid 'ncol'.", i);
            }
            if (TYPEOF(r_n) != INTSXP || LENGTH(r_n) != 1) {
                Rf_error("smats[%d] has invalid 'n'.", i);
            }
            if (TYPEOF(r_i) != INTSXP) {
                Rf_error("smats[%d] has invalid 'i'.", i);
            }
            if (TYPEOF(r_j) != INTSXP) {
                Rf_error("smats[%d] has invalid 'j'.", i);
            }
            if (TYPEOF(r_x) != REALSXP) {
                Rf_error("smats[%d] has invalid 'x'.", i);
            }

            data.smats[i] = Calloc(1, inla_cgeneric_smat_tp);
            data.smats[i]->nrow = INTEGER(r_nrow)[0];
            data.smats[i]->ncol = INTEGER(r_ncol)[0];
            data.smats[i]->n = INTEGER(r_n)[0];
            data.smats[i]->i = INTEGER(r_i);
            data.smats[i]->j = INTEGER(r_j);
            data.smats[i]->x = REAL(r_x);

            if (smat_names != R_NilValue && LENGTH(smat_names) > i) {
                data.smats[i]->name = strdup(CHAR(STRING_ELT(smat_names, i)));
            } else {
                data.smats[i]->name = strdup("unnamed");
            }

            Rprintf("smats[%d] (%s): nrow=%d, ncol=%d, n=%d\n", i, data.smats[i]->name, data.smats[i]->nrow, data.smats[i]->ncol, data.smats[i]->n);
            for (int j = 0; j < data.smats[i]->n && j < 10; j++) {
                Rprintf("  i[%d] = %d, j[%d] = %d, x[%d] = %f\n", j, data.smats[i]->i[j], j, data.smats[i]->j[j], j, data.smats[i]->x[j]);
            }
            if (data.smats[i]->n > 10) {
                Rprintf("  ... (only showing the first 10 entries)\n");
            }
        }
    } else {
        data.smats = NULL; // No sparse matrices present
    }

    Rprintf("All data successfully parsed. Calling the dynamic function now...\n");
    
     // Define integer constants for commands
    #define CMD_Q 1
    #define CMD_GRAPH 2
    #define CMD_MU 3
    #define CMD_INITIAL 4
    #define CMD_LOG_NORM_CONST 5
    #define CMD_LOG_PRIOR 6

    // Dynamically call the loaded function
    double *result = dynamic_function(cmd_parsed, theta, &data);

    int M; // Variable to store the number of non-zero entries or other relevant sizes
    SEXP r_result;

    switch (cmd_parsed) {
        case CMD_Q: // Allocate for M + 2
            M = (int) result[1]; // Second element of result contains M
            r_result = PROTECT(allocVector(REALSXP, M + 2));
            for (int i = 0; i < M + 2; i++) {
                REAL(r_result)[i] = result[i];
            }
            break;

        case CMD_GRAPH: // Allocate for 2*M + 2
            M = (int) result[1]; // Second element of result contains M
            r_result = PROTECT(allocVector(REALSXP, 2 * M + 2));
            for (int i = 0; i < 2 * M + 2; i++) {
                REAL(r_result)[i] = result[i];
            }
            break;

        case CMD_MU: // Allocate for 1
        case CMD_LOG_NORM_CONST:
        case CMD_LOG_PRIOR:
            r_result = PROTECT(allocVector(REALSXP, 1));
            REAL(r_result)[0] = result[0];
            break;

        case CMD_INITIAL: // Allocate for the length of theta
            M = (int) result[0];
            r_result = PROTECT(allocVector(REALSXP, M+1));
            for (int i = 0; i < M+1; i++) {
                REAL(r_result)[i] = result[i];
            }
            break;

        default: // Handle unknown command
            dlclose(handle);
            Rf_error("Unknown command received: %d", cmd_parsed);
    }

    // Release memory and return
    UNPROTECT(1);

    // Free allocated memory
    if (data.ints) {
        for (int i = 0; i < data.n_ints; i++) {
            Free(data.ints[i]->name);
            Free(data.ints[i]);
        }
        Free(data.ints);
    }
    if (data.doubles) {
        for (int i = 0; i < data.n_doubles; i++) {
            Free(data.doubles[i]->name);
            Free(data.doubles[i]);
        }
        Free(data.doubles);
    }
    if (data.chars) {
        for (int i = 0; i < data.n_chars; i++) {
            if (data.chars[i]) {
                Free(data.chars[i]->name);    // Free the name string
                Free(data.chars[i]->chars);  // Free the single string
                Free(data.chars[i]);         // Free the struct itself
            }
        }
        Free(data.chars);  // Free the array of struct pointers
    }
    if (data.mats) {
        for (int i = 0; i < data.n_mats; i++) {
            Free(data.mats[i]->name);
            Free(data.mats[i]);
        }
        Free(data.mats);
    }
    if (data.smats) {
        for (int i = 0; i < data.n_smats; i++) {
            Free(data.smats[i]->name);
            Free(data.smats[i]);
        }
        Free(data.smats);
    }

    dlclose(handle);
    return r_result;
}