#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_linalg.h>


int size;
// covariance matrix
double *covar;
// inversion of the covariance matrix
double *covar_inv;

int Ndim=2;

void likelihood_init()
{
	int i=0, j=0;	

	// allocate memory for the covariance matrix
	size = Ndim;
	covar     = malloc(sizeof(double)*size*size);
	covar_inv = malloc(sizeof(double)*size*size);

	// associate covar & covar_inv with gsl matrices
	gsl_matrix_view gsl_covar    = gsl_matrix_view_array(covar,     size, size);
	gsl_matrix_view gsl_covar_inv= gsl_matrix_view_array(covar_inv, size, size);
	
	// initialize covariance matrix
	for(i=0; i<size; i++)
	for(j=0; j<size; j++)
	{
		if(i != j)
			gsl_matrix_set(&gsl_covar.matrix, i, j, 0.5);
		else
			gsl_matrix_set(&gsl_covar.matrix, i, j, i+1.0);
	}

	// inverse the matrix
	int s;
	gsl_permutation * p = gsl_permutation_alloc(size);
	gsl_linalg_LU_decomp(&gsl_covar.matrix, p, &s);
	gsl_linalg_LU_invert(&gsl_covar.matrix, p, &gsl_covar_inv.matrix);

	// covar is not used any more
	free(covar);
}

double likelihood( int ndim, double *pos, int flag)
{
	double val=0.0;
	int i=0, j=0;

	for(i=0; i<ndim; i++)
	for(j=0; j<ndim; j++)
	{
		val += -0.5 * covar_inv[i*ndim+j]*pos[i]*pos[j];
	}
//	val = fmax(val, -MAXDOUBLE);
	return val;
}

void likelihood_finalize()
{
	free(covar_inv);
}
