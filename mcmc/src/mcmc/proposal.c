#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "proposal.h"

struct proposal_info Props;

/*=================================================================================================*/
void alloc_proposal_info(int ndim)
//	allocate proposal
{
	Props.ndim = ndim;
	Props.code =  malloc(sizeof(int)*ndim);
	Props.width = malloc(sizeof(double)*ndim);
}

void free_proposal_info()
{
	free(Props.code);
	free(Props.width);
}

void broadcast_proposal_info(void)
{
	int i;
	MPI_Aint indices[3];
	MPI_Datatype proposal_struct;

	int blocklens[3] = {1, Props.ndim, Props.ndim};
	MPI_Datatype old_types[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};

	MPI_Address( &Props.ndim,     &indices[0] );
	MPI_Address( &Props.code[0],  &indices[1] );
	MPI_Address( &Props.width[0], &indices[2] );

	int base = indices[0];
	for (i=0; i <3; i++) indices[i] -= base;

	//for (i=0; i<3; i++) printf("%ld \n", indices[i]);

	MPI_Type_struct( 3, blocklens, indices, old_types, &proposal_struct);
	MPI_Type_commit( &proposal_struct );

	MPI_Bcast( &Props, 1, proposal_struct, 0, MPI_COMM_WORLD );
}

/*========================================load proposal info from file=============================================*/

void read_proposal_info(char filename[])
{
	char buf[1000];
	int i, ibuf;
	FILE *fp;

	if( !(fp=fopen(filename, "r")))
	{
		printf("can't open file %s!\n", filename);
		exit(0);
	}

	fscanf(fp, "%d", &ibuf);
	if(ibuf == 0)
	{
		double a;
		int    c;
		fscanf(fp, "%d %lf %d", &ibuf, &a, &c);
		for(i=0; i<Props.ndim; i++)
		{
			Props.width[i] = a;
			Props.code[i]  = c;
		}
	}
	else
	{
		if(ibuf != Props.ndim )
		{
		        printf("inconsistent proposal file!\n");
			exit(1);
		}
		for(i=0;i<Props.ndim;i++)
		{
			fscanf(fp, "%d %lf %d", &ibuf, &Props.width[i], &Props.code[i]);
			if(Props.code[i] > 1)
			{
				printf("Wrong proposal function code: dim=%d code=%d", i, Props.code[i]);
				exit(1);
			}
		}
	}
	fclose(fp);
}

/*
void proposal_init(void)
//	called after run parameter struct is broadcasted
{
	int myid = MYID;

	alloc_proposal_info(Run.ndim);

	if(myid == 0)
		read_proposal_info(PropFileName);

	broadcast_proposal_info();	
}
*/

/*=========================================checkpoint==============================================*/
void fwrite_proposal_info(FILE *stream)
{
	fwrite(Props.code,  Props.ndim, sizeof(int),    stream);
	fwrite(Props.width, Props.ndim, sizeof(double), stream);
}

void fread_proposal_info(FILE *stream)
{
	fread(Props.code,  Props.ndim, sizeof(int),    stream);
	fread(Props.width, Props.ndim, sizeof(double), stream);
}
