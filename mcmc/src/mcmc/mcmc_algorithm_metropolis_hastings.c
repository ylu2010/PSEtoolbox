/*
 *  mcmc_algorithm_metropolis_hastings.c
 *  
 *
 *  Created by Y Lu on 4/12/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "run_info.h"
#include "proposal.h"
#include "chains.h"
#include "probability.h"
#include "mpi_rng.h"
#include "mcmc.h"

void metropolis_hastings_proposal(int ichain, int nvec, int ndim, double *state, double *t_state, double factor)
/*
	factor:	1 for MH
*/
{
	int idim;
	double dx;
	
	for(idim=0; idim<ndim; idim++)
	{
		dx = Props.width[idim] * random_number(Props.code[idim]) * factor;
		t_state[idim] = state[idim] + dx;
	}
}

void propose_state( int mcmc_method, int ichain, int nvec, int ndim, double *state, double *t_state, double factor)
{
	switch (mcmc_method)
	{
		case 0:
			metropolis_hastings_proposal(ichain, nvec, ndim, state, t_state, factor);
			break;
		case 1:
			differential_evolution_proposal(ichain, Run.nchain, nvec, ndim, state, t_state, factor);
			break;
		default:
			printf("propose_state: wrong mcmc_method!mcmc_method=%d \n", mcmc_method);
			exit(0);
	}
}

int accept_reject(double prob, int ichain, int nvec, int ndim, double *state, double *t_state)
{
	int i;
	double rand = random_number(0);
		
	if( exp(prob) > rand )
	{
		for(i=0; i<nvec; i++)
		{
			state[i] = t_state[i];
		}
		return 1;
	}
	else
		return 0;
}


int one_mcmc_step(int mcmc_method, int ichain, int nvec, int ndim, double *state, double beta, double factor)
{
	int accept = 0;	
	double prob_new, prob, accept_lnprob;
	double *t_state;
	
	prob = state[ndim];
	
	t_state = malloc(nvec * sizeof(double));
	propose_state(mcmc_method, ichain, nvec, ndim, state, t_state, factor);
	
	probability_evaluation(ichain, nvec, ndim, t_state);
	
	prob_new = t_state[ndim];
	
	accept_lnprob = (prob_new - prob) * beta;
	
	accept = accept_reject(accept_lnprob, ichain, nvec, ndim, state, t_state);
	
	free(t_state);

	return accept;
}
