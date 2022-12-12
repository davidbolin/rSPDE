#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"
#include "cgeneric_test_aux.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

double *inla_cgeneric_cppmodel(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{
	// this reimplement `inla.rgeneric.iid.model` using cgeneric

	double *ret = NULL, prec = (theta ? exp(theta[0]) : NAN), lprec = (theta ? theta[0] : NAN);

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(!strcasecmp(data->smats[0]->name, "sp_mat"));
  	inla_cgeneric_smat_tp *sp_mat = data->smats[0];

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the total length while N is the dimension

		int M = N;
		ret = Calloc(2 + 2 * N, double);
		assert(ret);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (int i = 0; i < M; i++) {
			ret[2 + i] = i;			       /* i */
			ret[2 + N + i] = i;		       /* j */
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		int M = N;
		ret = Calloc(2 + N, double);
		assert(ret);
		ret[0] = -1;				       /* REQUIRED! */
		ret[1] = M;				       /* number of (i <= j) */

		compute_Q(N, ret, sp_mat->x, sp_mat->i, sp_mat->j, sp_mat->n);

		for (int i = 0; i < M; i++) {
			ret[2 + i] = prec;
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters
		ret = Calloc(2, double);
		assert(ret);
		ret[0] = 1;
		ret[1] = 4.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.9189385332 + 0.5 * lprec);
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = -prec + lprec;			       // prec ~ gamma(1,1)
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}
