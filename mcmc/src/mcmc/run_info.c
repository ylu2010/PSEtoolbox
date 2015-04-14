#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "parameter.h"
#include "run_info.h"

struct run_info Run;
struct mcmc_info Mcmc;

void set_run_info(int numprocs)
{
// run
	if( !(get_parameter("NChain", INT, &Run.nchain)) )
	{
		printf("please specify the number of chains!\n");
		exit(0);
	}
	if( !(get_parameter("NDimension", INT, &Run.ndim)) )
	{
		printf("please specify the dimensionality!\n");
		exit(0);
	}
	Run.nproc = numprocs;
	Run.noffset = 3;
	Run.nvec = Run.noffset + Run.ndim;
	if( !(get_parameter("MaxIteration", INT, &Run.max_iter)) )
	{
		printf("please specify the max number of iterations!\n");
		exit(0);
	}
	if( !(get_parameter("OutputFileName", STRING, Run.chainlog)) )
	{
		strcpy(Run.chainlog, "woodpecker.dat");
	}

// mcmc
	if( !(get_parameter("MCMCMethod", INT, &Mcmc.method)) )
	{
		printf("please specify the mcmc algorithm!\n");
		exit(0);
	}
	// 0 for mh
	if(Mcmc.method == 0)	// mh
	{
		if( !(get_parameter("MaxTemperature", DOUBLE, &Mcmc.max_temp)) )
			Mcmc.max_temp = 1.0;
		if( Mcmc.max_temp > 1.0 )	// tempered algorithm
		{
			if( !(
					get_parameter("TempFreq", INT, &Mcmc.temp_freq) &&
				get_parameter("MinTempLevel", INT, &Mcmc.min_temp_level) &&
				get_parameter("MaxTempLevel", INT, &Mcmc.max_temp_level) &&
				get_parameter("IterPerTemp",  INT, &Mcmc.iter_per_temp)  &&
				get_parameter("ProposalScaleTempIndex", DOUBLE, &Mcmc.prop_scale_temp_index) 
			) )
			{
				printf("tempered transition is requested, but not all the relevant parameters are properly assigned.\n");
				printf("TempFreq\nMinTempLevel\nMaxTempLevel\nIterPerTemp\nProposalScaleTempIndex\n");
				exit(0);
			}
		}
	}
	else if(Mcmc.method == 1)	// differential evolution
	{
		if(Run.nchain <= 2)
		{
			printf("Differential Evolution requires the number of chains > 2!\n");
			exit(1);
		}

		get_parameter("DEFullJumpFreq", INT, &Mcmc.de_full_jump_freq);
		if( !(get_parameter("DEProposalScale", DOUBLE, &Mcmc.de_prop_scale)) )
			Mcmc.de_prop_scale = 2.38/sqrt(2.*Run.ndim);

		if( !(get_parameter("MaxTemperature", DOUBLE, &Mcmc.max_temp)) )
			Mcmc.max_temp = 1.0;
		if( Mcmc.max_temp > 1.0 )	// tempered algorithm
		{
			if( !(
				get_parameter("TempFreq", INT, &Mcmc.temp_freq) &&
				get_parameter("MinTempLevel", INT, &Mcmc.min_temp_level) &&
				get_parameter("MaxTempLevel", INT, &Mcmc.max_temp_level) &&
				get_parameter("IterPerTemp",  INT, &Mcmc.iter_per_temp)  &&
				get_parameter("ProposalScaleTempIndex", DOUBLE, &Mcmc.prop_scale_temp_index) 
			) )
			{
				printf("tempered transition is requested, but not all the relevant parameters are properly assigned.\n");
				printf("TempFreq\nMinTempLevel\nMaxTempLevel\nIterPerTemp\nProposalScaleTempIndex\n");
				exit(0);
			}
		}
	}
}

void broadcast_run_info( void )
{
	MPI_Bcast( &Mcmc, sizeof(struct mcmc_info), MPI_BYTE, 0, MPI_COMM_WORLD );
	MPI_Bcast( &Run,  sizeof(struct run_info),  MPI_BYTE, 0, MPI_COMM_WORLD );
}

void fwrite_run_info(FILE *stream)
{
	fwrite(&Run, 1, sizeof(struct run_info), stream);
	fwrite(&Mcmc,1, sizeof(struct mcmc_info),stream);
}

void fread_run_info(FILE *stream)
{
	fread(&Run, 1, sizeof(struct run_info), stream);
	fread(&Mcmc,1, sizeof(struct mcmc_info),stream);
}
