#ifndef chains_h
#define chains_h

extern struct state_info
{
	double prob, prior, likelihood;
	double *coord;
} S;

extern double *Chains, *State, *Coor, Prob, Prio, Like;

void alloc_chains( void );
void free_chains( void );

void init_chains_from_file(char filename[]);
void init_chains_from_prior();
void copy_state(int n, double *p_o, double *p_d);

void fwrite_chains(FILE *stream);
void fread_chains(FILE *stream);

#endif
