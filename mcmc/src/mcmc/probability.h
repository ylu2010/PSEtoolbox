#ifndef probability_h
#define probability_h

void probability_init(void);
void probability_evaluation( int ichain, int nvec, int ndim, double *state );
void probability_finalize(void);

#endif
