#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "parameter.h"
#include "checkpoint_io.h"
#include "chain_log_io.h"
#include "run_info.h"
#include "proposal.h"
#include "prior.h"
#include "chains.h"
#include "mcmc_log.h"
#include "mpi_rng.h"

void start_mcmc(int myid, int nproc, char * run_parameter_file)
{
	if(myid == 0)
		read_parameter(run_parameter_file);

	checkpoint_init(0, myid, nproc);
	
	if(myid == 0)
		set_run_info(nproc);
	broadcast_run_info();
	
	// proposal
	alloc_proposal_info(Run.ndim);
	if(myid == 0)
	{
		char proposal[50];
		if( !(get_parameter("ProposalFileName", STRING, proposal)) )
		{
			printf("start_mcmc: please specify the proposal file.\n");
			exit(0);
		}
		read_proposal_info(proposal);
	}
	broadcast_proposal_info();

	// prior
	alloc_prior_info(Run.ndim);
	if(myid == 0)
	{
		char prior[50];
		if( !(get_parameter("PriorFileName", STRING, prior)) )
		{
			printf("start_mcmc: please specify the prior file.\n");
			exit(0);
		}
		read_prior_info(prior);
	}
	broadcast_prior_info();

	// random number generator
	init_gsl_rng(myid);

	// chains
	alloc_chains();
	// initialize chains
	if(myid == 0)
	{
		int InitialStateMethod;
		char InitialStateFile[50];
		if( !(get_parameter("InitialStateMethod", INT, &InitialStateMethod)) )
		{
			printf("start_mcmc: please specify the method for initializing the chains.\n");
			exit(0);
		}
		if(InitialStateMethod == 1)
		{
			init_chains_from_prior();
		}
		else if(InitialStateMethod == 0)
		{
			if( !(get_parameter("InitialStateFile", STRING, InitialStateFile)) )
			{
				printf("start_mcmc: please specify the initial state file.\n");
				exit(0);
			}
			init_chains_from_file(InitialStateFile);
		}
	}

	// chain log
	if(myid == 0)
		new_chain_log();

	
	//mcmc diagnostic
	mcmc_log_init(myid, nproc);
	if(myid == 0)
		new_mcmc_log();
}

void resume_mcmc(int myid, int nproc, int restart, char * run_parameter_file)
{
	FILE * stream;

	if(myid == 0)
		read_parameter(run_parameter_file);

	checkpoint_init(restart, myid, nproc);
	
	// open checkpoint file to read
	stream = checkpoint_open_to_read(myid, nproc);
	checkpoint_read(stream, myid, nproc);
	checkpoint_close(stream, myid, nproc);

	// trim the log file
	if(myid == 0)
		trim_chain_log(Run.chainlog, Run.nchain*Run.iter+1);

	// something can be changed, while others can not
	// max # of interations can be changed
	if(myid == 0)
	{	
		if(get_parameter("MaxIteration", INT, &Run.max_iter))
			fputs("max number of interations changed.\n", stdout);
	}
	MPI_Bcast(&Run.max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// unfortunately, nothing else can be changed so far.
}

void end_mcmc()
{
	// free proposal
	free_proposal_info();
	// free prior
	free_prior_info();
	// free gsl random number generator
	free_gsl_rng();
	// free chains
	free_chains();
	// mcmc diagnostic
	mcmc_log_final();
}
