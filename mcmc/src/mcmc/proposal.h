#ifndef proposal_h
#define proposal_h

extern struct proposal_info
{
	int ndim;
	int *code;
	double *width;
} Props;

void alloc_proposal_info(int ndim);
void free_proposal_info();
void broadcast_proposal_info(void);

void read_proposal_info(char filename[]);

void fwrite_proposal_info(FILE *stream);
void fread_proposal_info(FILE *stream);

#endif
