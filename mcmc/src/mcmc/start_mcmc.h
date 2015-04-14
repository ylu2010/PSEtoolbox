#ifndef start_mcmc_h
#define start_mcmc_h

void start_mcmc(int myid, int nproc, char * run_parameter_file);
void resume_mcmc(int myid, int nproc, int restart, char * run_parameter_file);
void end_mcmc();

#endif
