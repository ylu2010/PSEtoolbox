#ifndef chain_log_io_h
#define chain_log_io_h

void write_chain_log(void);
void new_chain_log(void);
void trim_chain_log(char *filename, unsigned int nline);

#endif
