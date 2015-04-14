#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DOUBLE 1
#define FLOAT 2
#define INT 3
#define STRING 9

#define MAXPARAM 100

char tag[MAXPARAM][100];
char val[MAXPARAM][100];
int n_param;	

void read_parameter(char *fname);
int get_parameter(char *name, int type, void *value);

void read_parameter(char *fname)
{

        /* buf1: parameter name;
           buf2: parameter value;
           buf3: comments
        */
  	FILE *fd;
  	char buf[400], buf1[400], buf2[400], buf3[400];
	
  	printf("\nreading parameter file:\n");

	n_param = 0;
  	if((fd = fopen(fname, "r")))
  	{
		while( !feof(fd) && n_param < MAXPARAM )
		{
	  	        *buf = 0;
  		        /* load in one line ?*/
        		fgets(buf, 200, fd);
			/* check how many strings are there in this line. */
	      		if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
				continue;
          		/* jump over the comments marked by % */
          		if(buf1[0] == '%')
				continue;
			strcpy(tag[n_param], buf1);
			strcpy(val[n_param], buf2);
			printf("%s\t%s\n", tag[n_param], val[n_param]);
			n_param++;
      		}
      		fclose(fd);
  	}
  	else
  	{
      		printf("read_parameter: Parameter file %s not found.\n", fname);
    		exit(1);
  	}

	if(n_param == MAXPARAM)
	{
		printf("read_parameter: not all parameters are loaded, please check!");
		exit(0);
	}

}

int get_parameter(char *name, int type, void *value)
{
	int i, j;

	/* check the parameter list */
	for(i = 0, j = -1; i < n_param; i++)
	if(strcmp(name, tag[i]) == 0)
	{
		j = i;
		break;
	}
	/* load the parameter value */
	if(j >= 0)
	{
		switch (type)
        	{
			case DOUBLE:
        			*((double *) value) = atof(val[j]);
                  		break;
			case FLOAT:
				*((float *) value) = atof(val[j]);
				break;
                	case INT:
                  		*((int *) value) = atoi(val[j]);
                  		break;
                	case STRING:
                  		strcpy(value, val[j]);
                  		break;
		}
		return 1;
	}
	else
	{
		return 0;
        }
}
