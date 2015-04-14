#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "run_info.h"
#include "proposal.h"
#include "prior.h"
#include "mpi_rng.h"
#include "chains.h"
#include "mcmc_log.h"
#include "parameter.h"

static char checkpoint_base[50];
static unsigned int  checkpoint_code;
static unsigned int  checkpoint_interval;

void checkpoint_init(int restart, int myid, int nproc)
{
	if(myid == 0)
	{
		if( !(get_parameter("checkpoint_base",  STRING, checkpoint_base)) )
		{
			printf("checkpoint_init: please specify name base for checkpoint files.\n");
			exit(0);
		}
		if( !(get_parameter("checkpoint_interval", INT, &checkpoint_interval)) )
		{
			printf("checkpoint_init: checkpoint interval not specified. \n");
			printf("checkpoint_init: a large number is assigned.\n");
			checkpoint_interval = 4294967295;
		}
		checkpoint_code = restart;
	}
	MPI_Bcast(&checkpoint_interval, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int is_to_checkpoint(int iter, int myid, int nproc)
{
	if(iter%checkpoint_interval == 0)
		return 1;
	else
		return 0;
}

FILE * checkpoint_open_to_write(int myid, int nproc)
{
	char filename[70];

	if(myid == 0)
	{
		checkpoint_code += 1;
		sprintf(filename, "%s_chk_%.3d", checkpoint_base, checkpoint_code);
		return fopen(filename, "wb");
	}
	else
	{
		return NULL;
	}
}

FILE * checkpoint_open_to_read(int myid, int nproc)
{
	char filename[70];

	if(myid == 0)
	{
		sprintf(filename, "%s_chk_%.3d", checkpoint_base, checkpoint_code);
		return fopen(filename, "rb");
	}
	else
	{
		return NULL;
	}
}

void checkpoint_close(FILE *stream, int myid, int nproc)
{
	if(myid == 0)
		fclose(stream);
}

void checkpoint_write(FILE * stream, int myid, int nproc)
{

	if(myid == 0)
	{
		// run parameter
		fwrite_run_info(stream);
		// proposal
		fwrite_proposal_info(stream);
		// prior
		fwrite_prior_info(stream);
	}

	// random number generator
	checkpoint_gsl_rng(stream, 0, MPI_COMM_WORLD);

	if(myid == 0)
	{
		// chains
		fwrite_chains(stream);
		// mcmc diagnostic
		fwrite_mcmc_log(stream);
	}
}

void checkpoint_read(FILE *stream, int myid, int nproc)
{
	int nproc_prev;

	// run info
	if(myid == 0)
	{
		fread_run_info(stream);
	}
	broadcast_run_info();

	nproc_prev = Run.nproc;
	Run.nproc = nproc;

	// proposal
	alloc_proposal_info(Run.ndim);
	if(myid == 0)
	{
		fread_proposal_info(stream);
	}
	broadcast_proposal_info();

	// prior
	alloc_prior_info(Run.ndim);
	if(myid == 0)
	{
		fread_prior_info(stream);
	}
	broadcast_prior_info();

	// random number generator
	// allocate memory first!
	init_gsl_rng(myid);
	resume_gsl_rng(stream, nproc_prev, 0, MPI_COMM_WORLD);

	// load chains
	alloc_chains();
	if(myid == 0)
	{
		fread_chains(stream);
	}

	// mcmc diagnostic
	mcmc_log_init(myid, nproc);
	if(myid == 0)
		fread_mcmc_log(stream);

}
