/*
 *  mcmc_algorithm_tempered_transition.c
 *  
 *	Implementation of the Tempered Transitions by Radford Neal 
 *	Neal. Sampling from multimodal distributions using tempered transitions. Stat Comput (1996)
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
#include "mpi_rng.h"
#include "mcmc.h"

struct tempred_transition_setting
{
	int ntemp;
	int ninter;
	double *beta;
	double *factor;
	int cur_level;
} tts;

double getTemperature()
{
	if(Mcmc.max_temp > 1.) {
	    return 1./tts.beta[tts.cur_level];
	} else {
	    return 1.;
	}
}
	
void tempered_transition_init(int ninter, int minmc, int maxmc, int ndim, int itmax, double tmax, double tpow)
// ninter:number of interation on each temperature
// minmc: mimimum number of temperature levels
// maxmc: maximum number of temperature levels
// ndim:  dimension of the problem
// itmax:
// tpow:  0.5 for Gaussian scaling
// tmax:  maximum temperature
{
	int k;
	
	tts.ntemp  = fmax( minmc, (int) (sqrt(ndim + 3.0) * log(tmax)));
	tts.ntemp  = fmin( maxmc, tts.ntemp);
	tts.ninter = ninter;
	tts.beta   = malloc((tts.ntemp+1) * sizeof(double));
	tts.factor = malloc((tts.ntemp+1) * sizeof(double));
	
	printf("Initializing Tempered Transitions with %d temperature levels and %d iteration per temperature.\n", tts.ntemp, tts.ninter);
	tts.beta[0] = tts.factor[0] = 1.;
	for(k=1; k<=tts.ntemp; k++)
	{
		tts.beta[k]=exp(-log(tmax)*k/(tts.ntemp));
		tts.factor[k] = 1./pow(tts.beta[k], tpow);
	}
	tts.cur_level = 0;
}

int tempered_transition_step(int method, int ichain, int nvec, int ndim, double *state)
{
	int i, j, k;
	double up=0.0, down=0.0;
	double *t_state;
	double en, accept_lnprob;
	int accept = 0;
	
	t_state = malloc(nvec*sizeof(double));
	for (j=0; j<nvec; j++) t_state[j]=state[j];
	
	for (k=1; k<=tts.ntemp; k++) 
	{
		tts.cur_level = k;
		en = t_state[ndim];
		up += en * (tts.beta[k] - tts.beta[k-1]);
		for (i=0; i<tts.ninter; i++) one_mcmc_step(method, ichain, nvec, ndim, t_state, tts.beta[k], tts.factor[k]); //mcmc_method =0
	}
	
	for (k=tts.ntemp; k>=1; k--)
	{
		tts.cur_level = k;
		for (i=0; i<tts.ninter; i++) one_mcmc_step(method, ichain, nvec, ndim, t_state, tts.beta[k], tts.factor[k]);
		en = t_state[ndim];
		down += en * (tts.beta[k] - tts.beta[k-1]);
	}
	// reset the current temperature level
	tts.cur_level = 0;
	
	accept_lnprob = fmin(0.0, (up - down));
	accept = accept_reject(accept_lnprob, ichain, nvec, ndim, state, t_state);
	
	free(t_state);
	
	return accept;
}

void tempered_transition_end(void)
{
	free(tts.beta);
	free(tts.factor);
}

