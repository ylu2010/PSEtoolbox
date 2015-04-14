#ifndef parameter_h
#define parameter_h

#define DOUBLE 1
#define FLOAT 2
#define INT 3
#define STRING 9

void read_parameter(char *fname);
int get_parameter(char *name, int type, void *value);

#endif
