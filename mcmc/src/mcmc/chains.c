#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "run_info.h"
#include "prior.h"
#include "chains.h"

struct state_info S;

double *Chains, *State, *Coor, Prob, Prio, Like;

/*=================================================================================================*/
void alloc_chains( void )
//	allocate chain states
{
	Chains = malloc(sizeof(double) * Run.nvec * Run.nchain);
	State =  malloc(sizeof(double) * Run.nvec);
}

void free_chains(void)
{
	free(Chains);
	free(State);
}

void init_chains_from_file(char filename[])
{
	int i, j, k, indx, ibuf;
	FILE *fp;

	// start from initial file states
	if( !(fp=fopen(filename, "r")))
	{
		printf("chains: can't open file %s!\n", filename);
		exit(0);
	}

	fscanf(fp, "%d", &ibuf);
	if(ibuf != Run.ndim )
	{
		printf("chains: inconsistent initial state file!\n");
		exit(0);
	}
	fscanf(fp, "%d", &ibuf);
	if(ibuf != Run.nchain )
	{
		printf("chains: inconsistent initial state file!\n");
		exit(0);
	}
	for(i=0; i<Run.nchain; i++)
	{
		k = i*Run.nvec;
		fscanf(fp, "%d ", &ibuf);
		for(j=0, indx=0;j<Run.ndim;j++)
		{
      			fscanf(fp, "%lf", &Chains[k+j]);
		}
	}
	fclose(fp);
}

void init_chains_from_prior()
{
	int i, j, k;
	double v;

	for(i=0; i<Run.nchain; i++)
	{
		k = i*Run.nvec;
		for(j=0; j<Run.ndim; j++)
		{
			v = random_sample_prior( j );
			Chains[k+j] = v;
		}
	}
}

void copy_state(int n, double *p_o, double *p_d)
{
	int i;
	for(i=0; i<n; i++) p_d[i] = p_o[i];
}

/*=============================================checkpoint================================================*/
void fwrite_chains(FILE *stream)
{
	fwrite(Chains, Run.nvec*Run.nchain, sizeof(double), stream);
}

void fread_chains(FILE *stream)
{
	fread(Chains, Run.nvec*Run.nchain, sizeof(double), stream);
}
