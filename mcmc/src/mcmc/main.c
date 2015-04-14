#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#include "start_mcmc.h"
#include "checkpoint_io.h"
#include "mpi_rng.h"
#include "mcmc.h"
#include "probability.h"

int main(int argc, char *argv[])
{
	int restart = 0;
	char run_parameter_file[100] = "PARAMETER";

	int myid, numprocs;
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(processor_name, &namelen);

	while(1)
	{
		int c = getopt(argc, argv, "r:f:");
		if( c == -1 ) break;
		switch(c)
		{
			case 'r': restart = atoi(optarg);	 	break;
			case 'f': strcpy(run_parameter_file, optarg);	break;
			case '?': 
				fputs("please offer values for each of the options.\n", stderr);
				abort();
		}	
	}

	if(restart == 0)	// start a new mcmc simulation
	{
		printf("starting a new simulation...\n");
		start_mcmc(myid, numprocs, run_parameter_file);
	}
	else			// resume a previous simulation
	{
		printf("resuming a simulation from checkpoint file #%d...\n", restart);
		resume_mcmc(myid, numprocs, restart, run_parameter_file);
	}
	printf("\n");

	probability_init();

	run_mcmc(myid, numprocs);

	end_mcmc();

	MPI_Finalize();
	return 1;
}
