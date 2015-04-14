#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "run_info.h"
#include "chains.h"
#include "proposal.h"
#include "chain_log_io.h"
#include "checkpoint_io.h"
#include "mpi_rng.h"

#include "mcmc.h"
#include "mcmc_log.h"

void run_mcmc(int myid, int nproc)
{

	int nloop, iloop, ichain, p, p0;
	double *s_temp;
	double *state;
	int accept;

	if (Mcmc.max_temp > 1.)
	{
		printf("Tempered transitions is requested, Tmax = %g\n", Mcmc.max_temp);
		tempered_transition_init( Mcmc.iter_per_temp, Mcmc.min_temp_level, Mcmc.max_temp_level, Run.ndim, 100, Mcmc.max_temp, Mcmc.prop_scale_temp_index);
	}
	
	state = malloc(Run.nvec * sizeof(double));
	
	nloop = (int) (Run.nchain/Run.nproc);

	MPI_Bcast( Chains, Run.nvec*Run.nchain, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	// initialize the probabilities
	for(iloop=0; iloop<nloop; iloop++)
	{
		ichain = Run.nproc * iloop + myid;
		p0 = Run.nproc * iloop * Run.nvec;
		p = ichain * Run.nvec;

		MPI_Scatter( &Chains[p0], Run.nvec, MPI_DOUBLE, state, Run.nvec, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		probability_evaluation( ichain, Run.nvec, Run.ndim, state );

		MPI_Gather( state, Run.nvec, MPI_DOUBLE, &Chains[p], Run.nvec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	// chain log 
	if(myid == 0)
	{
		write_chain_log();
	}

	Run.stop = 0;
	while(Run.stop == 0)
	{
		Run.iter++;
		random_number_uniform(0.0,1.0);

		MPI_Bcast( Chains, Run.nchain * Run.nvec, MPI_DOUBLE, 0, MPI_COMM_WORLD );

        	for(iloop=0; iloop<nloop; iloop++)
		{
			ichain = Run.nproc * iloop + myid;
			p = ichain * Run.nvec;
//
//			propose_state( Mcmc_method, ichain, Run.nvec, Run.ndim, state, t_state, 1 );
//			accept_state( ichain );
//
//			s_temp = &Chains[p];
//			MPI_Gather( s_temp, Run.nvec, MPI_DOUBLE, &Chains[p], Run.nvec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			copy_state(Run.nvec, &Chains[p], state);
			
			if(Mcmc.max_temp > 1. && !(Run.iter % Mcmc.temp_freq))
				accept = tempered_transition_step(Mcmc.method, ichain, Run.nvec, Run.ndim, state);
			else
				accept = one_mcmc_step( Mcmc.method, ichain, Run.nvec, Run.ndim, state, 1, 1); //beta = 1, factor = 1
			
			copy_state(Run.nvec, state, &Chains[p]);
		
			MPI_Gather( state, Run.nvec, MPI_DOUBLE, &Chains[p], Run.nvec, MPI_DOUBLE, 0, MPI_COMM_WORLD );

			// mcmc diagnositc
			MPI_Gather( &accept, 1, MPI_INT, &accept_1step[Run.nproc*iloop + myid], 1, MPI_INT, 0, MPI_COMM_WORLD );
		}
		// mcmc diagnositc
		if(myid == 0)
		{
			int i=0;
			for(i=0; i<Run.nchain; i++)
				accept_count[i] += accept_1step[i];
			write_mcmc_log();
		}

		// chain log 
		if(myid == 0)
		{
			write_chain_log();
		}

		// checkpoint
		if(is_to_checkpoint(Run.iter, myid, nproc))
		{
			FILE *stream;
			stream = checkpoint_open_to_write(myid, nproc);
			checkpoint_write(stream, myid, nproc);
			checkpoint_close(stream, myid, nproc);
		}

		if (Run.iter >= Run.max_iter) Run.stop = 1;
	}
	
	free(state);

	if(Mcmc.max_temp > 1.)
		tempered_transition_end();
}
