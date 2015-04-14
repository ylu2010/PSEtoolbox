#ifndef likelihood_h
#define likelihood_h

void likelihood_init();
double likelihood(int ndim, double * pos, int flag);
void likelihood_finalize();

#endif
