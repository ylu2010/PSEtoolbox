#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "run_info.h"
#include "proposal.h"
#include "chains.h"
#include "mpi_rng.h"
#include "mcmc.h"

void differential_evolution_proposal(int ichain, int nchain, int nvec, int ndim, double *state, double *t_state, double factor)
/*
	factor:	1 for DE
*/
{
	int ich_r1, ich_r2, ielem_r1, ielem_r2;
	int idim, ielem; 
	double gamma, dx;
	double *s_r1, *s_r2;

	int i = Run.iter % Mcmc.de_full_jump_freq;
	if( i ) 
		gamma = Mcmc.de_prop_scale * factor;
	else 
		gamma = 1.0;

	// ich_r1 ich_r2 & ichain are different from each other
	ich_r1 = floor(random_number(0)*(nchain-1));
	if(ich_r1 >= ichain) ich_r1++;
	ich_r2 = floor(random_number(0)*(nchain-2));
	if(ich_r2 >= ichain) ich_r2++;
	if(ich_r2 >= ich_r1) ich_r2++;

	s_r1 = &Chains[ich_r1*nvec];
	s_r2 = &Chains[ich_r2*nvec];
	for(idim = 0; idim < ndim; idim++)
	{
        	dx = gamma * (s_r2[idim] - s_r1[idim]) + Props.width[idim]*2.*(random_number(0)-0.5) * factor;
		t_state[idim] = state[idim] + dx;
	}
}
