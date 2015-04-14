#ifndef prior_h
#define prior_h

struct prior_info
{
	int ndim;
	int *code;
	double *min, *max;
};

extern struct prior_info Prior;

void alloc_prior_info(int ndim);
void free_prior_info();
void broadcast_prior_info(void);

void read_prior_info(char filename[]);

void fwrite_prior_info(FILE *stream);
void fread_prior_info(FILE *stream);

double random_sample_prior( int idim );
double distribution_uniform(int ndim, double *pos, struct prior_info prior);
double prior( int ndim, double *pos );
double ln_prior( int ndim, double *pos );

#endif
