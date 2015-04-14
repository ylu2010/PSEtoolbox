#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "likelihood.h"
#include "prior.h"

void probability_init(void)
{
	likelihood_init();
}

void probability_evaluation( int ichain, int nvec, int ndim, double *state )
{
	state[ndim+1] = ln_prior( ndim, state );
	state[ndim+2] = likelihood( ndim, state, ichain);
	state[ndim]   = state[ndim+1] + state[ndim+2];
}

void probability_finalize(void)
{
	likelihood_finalize();
}
