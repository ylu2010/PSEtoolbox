#ifndef run_info_h
#define run_info_h

#define MINDOUBLE 1e-33
#define MAXDOUBLE 1e33

extern struct run_info
{
	int nchain;
	int ndim;
	int nproc;
	int noffset;
	int nvec;
	int max_iter;
	char chainlog[50];
	int stop;
	int iter;
} Run;

extern struct mcmc_info
{
	int method;
	int temp_freq; // for tempered simulation
	int min_temp_level;
	int max_temp_level;
	int iter_per_temp;
	int de_full_jump_freq;		// for differential evolution
	double max_temp;
	double prop_scale_temp_index;
	double de_prop_scale;		// for differential evolution
} Mcmc;

void set_run_info(int numprocs);
void broadcast_run_info( void );

void fwrite_run_info(FILE *stream);
void fread_run_info(FILE *stream);

#endif
