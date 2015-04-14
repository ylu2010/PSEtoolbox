#ifndef mcmc_log_h
#define mcmc_log_h

extern unsigned int * accept_count;
extern unsigned int * accept_1step;

void mcmc_log_init(int myid, int nproc);
void mcmc_log_final(void);
void write_mcmc_log(void);
void new_mcmc_log(void);
void fwrite_mcmc_log(FILE *stream);
void fread_mcmc_log(FILE *stream);

#endif
