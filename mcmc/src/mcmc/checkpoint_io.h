#ifndef checkpoint_io_h
#define checkpoint_io_h

void checkpoint_init(int restart, int myid, int nproc);

int is_to_checkpoint(int iter, int myid, int nproc);

FILE * checkpoint_open_to_write(int myid, int nproc);
FILE * checkpoint_open_to_read(int myid, int nproc);
void checkpoint_close(FILE *stream, int myid, int nproc);

void checkpoint_write(FILE * stream, int myid, int nproc);
void checkpoint_read(FILE *stream, int myid, int nproc);

#endif
