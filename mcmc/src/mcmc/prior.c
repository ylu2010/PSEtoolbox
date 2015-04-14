#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "prior.h"
#include "mpi_rng.h"

#define MAXDOUBLE	1.e+308
#define MINDOUBLE	1.e-308

struct prior_info Prior;

//==================================================================================================

void alloc_prior_info(int ndim)
//	allocate prior
{
	Prior.ndim = ndim; 
	Prior.code = malloc(sizeof(int)*ndim);
	Prior.min  =  malloc(sizeof(double)*ndim);
	Prior.max  =  malloc(sizeof(double)*ndim);
}

void free_prior_info()
{
	free(Prior.code);
	free(Prior.min);
	free(Prior.max);
}

void broadcast_prior_info(void)
{
	int i;
	MPI_Aint indices[4];
	MPI_Datatype prior_struct;

	int blocklens[4] = {1, Prior.ndim, Prior.ndim, Prior.ndim};
	MPI_Datatype old_types[4] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Address( &Prior.ndim,    &indices[0] );
	MPI_Address( &Prior.code[0], &indices[1] );
	MPI_Address( &Prior.min[0],  &indices[2] );
	MPI_Address( &Prior.max[0],  &indices[3] );

	int base = indices[0];
	for (i=0; i <4; i++) indices[i] -= base;

	//for (i=0; i<4; i++) printf("%ld \n", indices[i]);

	MPI_Type_struct( 4, blocklens, indices, old_types, &prior_struct);
	MPI_Type_commit( &prior_struct );

	MPI_Bcast( &Prior, 1, prior_struct, 0, MPI_COMM_WORLD );
}

/*===============================initialize prior====================================*/

void read_prior_info(char filename[])
//	allocate prior_info first!
{
	int i, ibuf;
	FILE *fp;

	if( !(fp=fopen(filename, "r")))
	{
		printf("prior: can't open file %s!\n", filename);
		exit(0);
	}

	fscanf(fp, "%d", &ibuf);
	if(ibuf == 0)
	{
		double a, b;
		int    c;
		fscanf(fp, "%d %lf %lf %d", &ibuf, &a, &b, &c);
		for(i=0; i<Prior.ndim; i++)
		{
			Prior.min[i] = a;
			Prior.max[i] = b;
			Prior.code[i]= c;
		}
	}
	else
	{
		if(ibuf != Prior.ndim )
		{
			printf("prior: inconsistent prior file!\n");
			exit(1);
		}
		for(i=0;i<Prior.ndim;i++)
		{	
			fscanf(fp, "%d %lf %lf %d", &ibuf, &Prior.min[i], &Prior.max[i], &Prior.code[i]);
		}
	}
	fclose(fp);
}

/*
void prior_init(int myid)
//	called after run parameter struct is broadcasted
{
	alloc_prior_info(Run.ndim);

	if(myid == 0)
		read_prior_info(PriorFileName);

	broadcast_prior_info();	
}
*/

/*==============================checkpoint===============================*/

void fwrite_prior_info(FILE *stream)
{
	fwrite(Prior.code, Prior.ndim, sizeof(int),    stream);
	fwrite(Prior.min,  Prior.ndim, sizeof(double), stream);
	fwrite(Prior.max,  Prior.ndim, sizeof(double), stream);
}

void fread_prior_info(FILE *stream)
//	allocate prior_info first!
{
	fread(Prior.code, Prior.ndim, sizeof(int),    stream);
	fread(Prior.min,  Prior.ndim, sizeof(double), stream);
	fread(Prior.max,  Prior.ndim, sizeof(double), stream);
}

/*
void resume_prior_info(FILE *stream, int myid) 
//	called after run parameter struct is broadcasted
{
	alloc_prior_info(Run.ndim);	

	if(myid == 0)
		fread_prior_info(stream);

	broadcast_prior_info();	
}
*/

/*==============================sample from prior====================================*/

double random_sample_prior( int idim )
{
	double v;

	switch(Prior.code[idim])
	{
		case 0:
			v = random_number_uniform( Prior.min[idim], Prior.max[idim]);
			break;
		case 1:
			printf("gaussian prior has not been implemented!\n");
			exit(0);
			break;
		default:
			printf("wrong prior code %d!\n", Prior.code[idim]);
			exit(0);
	}
	return v;
}

/*===============================evaluate prior========================================*/

double distribution_uniform(int ndim, double *pos, struct prior_info prior)
{
	int i;
	double vol, l, p;

	for(i=0, vol=1.0; i<ndim; i++)
	{
		l = prior.max[i] - prior.min[i];
		vol *= l;
	}
	p = 1./vol;

	for(i=0, l=1; i<ndim; i++)
	{
		if(pos[i] < prior.min[i] || pos[i] > prior.max[i])
			l *= 0;
		else 
			l *= 1.;
	}
	if(l == 0) p = 0;

	return p;
}

double distribution_gaussian(int ndim, double *pos){return 0.0;}

double prior( int ndim, double *pos )
{
	double p;
	
	p = distribution_uniform( ndim, pos, Prior );

	return p;
}


double ln_prior( int ndim, double *pos )
{
	double p, lnp;

	p = prior( ndim, pos );
	if (p > MINDOUBLE) 
		lnp = log(p);
	else 
		lnp = -MAXDOUBLE;

	return lnp;
}
