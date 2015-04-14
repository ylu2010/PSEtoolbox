#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "run_info.h"
#include "chains.h"

void keep_record_if_max( int ich )
{
    static double max_prob = -MAXDOUBLE;
    int j;
    FILE *fm;

    if(!(fm = fopen("maxlog", "a+")))
    {
        printf("maxlog: can't open maxlog!\n");
        exit(1);
    }

    if( Prob > max_prob)
    {
        max_prob = Prob;
        fprintf(fm, "%8d %8d", ich, Run.iter);
        fprintf(fm, " %16.8g %16.8g %16.8g", Prob, Prio, Like);
        for(j=0; j<Run.ndim; j++)
        {
            fprintf(fm, " %16.8g", Coor[j]);
        }
        fprintf(fm, "\n");
    }
    fclose(fm);
}

void convert_chain_states( int p )
{
	Coor = &Chains[p];
	Prob = Chains[p+Run.ndim];
	Prio = Chains[p+Run.ndim+1];
	Like = Chains[p+Run.ndim+2];
}

void write_chain_log(void)
{
	char buf[50];
	int i, j, p;
	double lnpro, lnlik, lnpri;
	FILE *fd;

	sprintf(buf, "%s", Run.chainlog);

	if(!(fd = fopen(buf, "a+")))
	{
		printf("chain_log: can't open output file `%s' to write!\n", buf);
		exit(1);
	}

	for(i=0; i<Run.nchain; i++)
	{
		p = i*Run.nvec;
		convert_chain_states( p );
		keep_record_if_max( i );
		fprintf(fd, "%8d %8d", 0, Run.iter);

		lnpro = Prob;
		lnpri = Prio;
		lnlik = Like;

		fprintf(fd, " %16.8g %16.8g %16.8g", lnpro, lnpri, lnlik);
		for(j=0; j<Run.ndim; j++)
		{
			fprintf(fd, " %16.8g", Coor[j]);
		}
		fprintf(fd, "\n");
	}

	fclose(fd);
}

void new_chain_log(void)
{
	char buf[50];
	int i;
	FILE *fd;

	sprintf(buf, "%s", Run.chainlog);

	if(!(fd = fopen(buf, "w")))
	{
		printf("chain_log: can't open output file `%s' to write!\n", buf);
		exit(1);
	}
	
	// header of the log file
	fprintf(fd, "%8s %8s", "Level", "Iter");
	fprintf(fd, " %16s %16s %16s", "Probability", "Prior", "Likelihood");
	for(i=0; i<Run.ndim; i++)
	{
		fprintf(fd, "         Param%03d", i);
	}
	fprintf(fd, "\n");

	fclose(fd);
}

void trim_chain_log(char *filename, unsigned int nline)
{
	char backupfile[50];
	char line[2048];
	unsigned int i=0;

	FILE * backup, * origin;

	sprintf(backupfile, "%s.bak", filename);
	if( !(backup = fopen(backupfile, "w")) )
	{
		printf("chain_log: can't open file `%s' to write!\n", backupfile);
		exit(1);
	}
	if( !(origin = fopen(filename, "r")) )
	{
		printf("chain_log: can't open file `%s' to read!\n", filename);
		exit(1);
	}
		
	while(fgets(line, 1024, origin))
		fputs(line, backup);
	fclose(backup);
	fclose(origin);

	backup = fopen(backupfile, "r");
	origin = fopen(filename, "w");
	for(i=0; i<nline; i++)
	{
		fgets(line, 1024, backup);
		fputs(line, origin);
	}
	fclose(backup);
	fclose(origin);
}
