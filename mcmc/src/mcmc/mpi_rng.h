#ifndef mpi_rng_h
#define mpi_rng_h

#include <gsl/gsl_rng.h>
#include <mpi.h>

void init_gsl_rng( int iseed );
void free_gsl_rng(void);

void checkpoint_gsl_rng(FILE *stream, int root, MPI_Comm comm);
void resume_gsl_rng(FILE *stream, int n_rng, int root, MPI_Comm comm);

int broadcast_gsl_rng(gsl_rng *r, int root, MPI_Comm comm);
int send_gsl_rng(gsl_rng *r, int dest, int tag, MPI_Comm comm);
int recv_gsl_rng(gsl_rng *r, int source, int tag, MPI_Comm comm, MPI_Status *status);

double random_number(int flag);
double random_number_uniform(double x1, double x2);

#endif
