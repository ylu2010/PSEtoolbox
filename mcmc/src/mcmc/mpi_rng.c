#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "mpi_rng.h"

gsl_rng *rng;

/*==========================start new rng================================*/
void init_gsl_rng( int iseed )
{
	gsl_rng_env_setup();

	const gsl_rng_type *T;

	T = gsl_rng_default;
	rng = gsl_rng_alloc(T);

	gsl_rng_set(rng, (unsigned long int)iseed+1);
}

void free_gsl_rng(void)
{
	gsl_rng_free(rng);
}

/*==============================checkpoint================================*/
void checkpoint_gsl_rng(FILE *stream, int root, MPI_Comm comm)
/*
	gather all the generators and put them into the checkpoint file
*/
{
	int pid, np;
	MPI_Comm_rank(comm, &pid);
	MPI_Comm_size(comm, &np);

	// a temperary generator
	gsl_rng *t;
	t = gsl_rng_alloc(gsl_rng_default);
	
	if(pid == root)
	{
		int i=0;
		for(i=0; i<np; i++)
		{
			if( i == root )
			{
				gsl_rng_fwrite(stream, rng);
			}
			else
			{
				MPI_Status status;
				recv_gsl_rng(t, i, 0, comm, &status);
				gsl_rng_fwrite(stream, t);
			}
		}
	}
	else
	{
		send_gsl_rng(rng, root, 0, comm);
	}
}

void resume_gsl_rng(FILE *stream, int n_rng, int root, MPI_Comm comm)
/*
DESCRIPTION:
	resume random number generator from checkpoint file
	generators must be initialized first prior to this function
ARGUMENTS:
	n_rng:	number of generators in the checkpoint file
*/
{
	int pid, np, np_prev, np_min;
	MPI_Comm_rank(comm, &pid);
	MPI_Comm_size(comm, &np);
	np_prev = n_rng;

	if(np >= np_prev) 
		np_min = np_prev;
	else
		np_min = np;

	// a temperary generator
	gsl_rng *t;
	t = gsl_rng_alloc(gsl_rng_default);

	if(pid == root)
	{
		int i=0;
		for(i=0; i<np_min; i++)
		{
			if(i == root)
			{
				gsl_rng_fread(stream, rng);
			}
			else
			{
				gsl_rng_fread(stream, t);
				send_gsl_rng(t, i, 0, comm);
			}
		}
		for(i=np; i<np_prev; i++)	// some generators are left in the checkpoint
		{
			gsl_rng_fread(stream, t);
		}
	}
	else if(pid < np_min)
	{
		MPI_Status status;
		recv_gsl_rng(rng, root, 0, comm, &status);
	}
}

/*==========================parallel operations===========================*/

int broadcast_gsl_rng(gsl_rng *r, int root, MPI_Comm comm)
{
	return MPI_Bcast(r->state, (int)(r->type->size), MPI_CHAR, root, comm);
}

int send_gsl_rng(gsl_rng *r, int dest, int tag, MPI_Comm comm)
{
	return MPI_Send(r->state, (int)(r->type->size), MPI_CHAR, dest, tag, comm);
}

int recv_gsl_rng(gsl_rng *r, int source, int tag, MPI_Comm comm, MPI_Status *status)
{
	return MPI_Recv(r->state, (int)(r->type->size), MPI_CHAR, source, tag, comm, status);
}

/*====================sample from different distributions=======================*/
double random_number(int flag)
{
	double rn;

	if(flag == 0) rn = gsl_rng_uniform(rng);

	if(flag == 1) rn = gsl_ran_ugaussian(rng);
	//printf("flag=%d rn=%g\n", flag, rn);
	return rn;
}

double random_number_uniform(double x1, double x2)
{
	double r = (x2-x1)*random_number(0)+x1;
	return r;
}
