#include <stdio.h>
#include <stdlib.h>

#include "run_info.h"

unsigned int * accept_count;
unsigned int * accept_1step;

void mcmc_log_init(int myid, int nproc)
{
	int i=0;
	
	accept_count = malloc(sizeof(unsigned int)*Run.nchain);
	for(i=0; i<Run.nchain; i++)
		accept_count[i] = 0;

	accept_1step = malloc(sizeof(unsigned int)*Run.nchain);
}

void mcmc_log_final(void)
{
	free(accept_count);
	free(accept_1step);
}

// io
void write_mcmc_log(void)
{
	int i = 0;
	FILE *fd;

	if(!(fd = fopen("woodpecker.log", "a+")))
	{
		printf("mcmc_log: can't open log file `%s' to write!\n", "woodpecker.log");
		exit(1);
	}

	for(i=0; i<Run.nchain; i++)
	{
		fprintf(fd, "%8d %8d", Run.iter, i);
		fprintf(fd, " %8d %8d %8.2g", accept_1step[i], accept_count[i], accept_count[i]/(double)Run.iter);
		fprintf(fd, "\n");
	}

	fclose(fd);
}

void new_mcmc_log(void)
{
	int i;
	FILE *fd;

	if(!(fd = fopen("woodpecker.log", "w")))
	{
		printf("mcmc_log: can't open log file `%s' to write!\n", "woodpecker.log");
		exit(1);
	}
	
	// header of the log file
	fprintf(fd, "%8s %8s", "Iter", "Chain");
	fprintf(fd, " %8s %8s %8s", "1Step", "Count", "Rate");
	fprintf(fd, "\n");

	fclose(fd);
}

//checkpoint

void fwrite_mcmc_log(FILE *stream)
{
	fwrite(accept_count, Run.nchain, sizeof(unsigned int), stream);
}

void fread_mcmc_log(FILE *stream)
{
	fread(accept_count, Run.nchain, sizeof(unsigned int), stream);
}
