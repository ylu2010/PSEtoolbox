#ifndef mcmc_h
#define mcmc_h

void run_mcmc(int myid, int nproc);

// mcmc_algorithm_metropolis_hastings.c
void metropolis_hastings_proposal(int ichain, int nvec, int ndim, double *state, double *t_state, double factor);
int accept_reject(double prob, int ichain, int nvec, int ndim, double *state, double *t_state);
int one_mcmc_step(int mcmc_method, int ichain, int nvec, int ndim, double *state, double beta, double factor);

// mcmc_algorithm_differential_evolution.c
void differential_evolution_proposal(int ichain, int nchain, int nvec, int ndim, double *state, double *t_state, double factor);

// mcmc_algorithm_tempered_transition.c
double getTemperature();
void tempered_transition_init(int ninter, int minmc, int maxmc, int ndim, int itmax, double tmax, double tpow);
int tempered_transition_step(int mcmc_method, int ichain, int nvec, int ndim, double *state);
void tempered_transition_end(void);

#endif
