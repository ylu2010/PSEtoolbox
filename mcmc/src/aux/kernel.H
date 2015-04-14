#ifndef kernel_h
#define kernel_h

#ifdef __cplusplus

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

using namespace std;

class Kernel
{
public:
	unsigned int ndim;
	// kernel center
	vector<double> coord;
	// pdf
	double pdf;
	// kernel size
	double size;

	// Gaussian kernel
	double gaussian(vector<double> x);
};

class Kernels
{
protected:
	vector<Kernel> kernel;
	unsigned int num_points;
	unsigned int ndim;

public:
	// constructor
	Kernels(const char * fname);

	// access
	unsigned int get_ndim()
	{
		return ndim;
	}
	unsigned int get_num()
	{
		return num_points;
	}

	// return pdf
	double Density(vector<double> x);
};

#else
typedef 
    struct Kernels
    Kernels;
#endif

#ifdef __cplusplus
extern "C" {
#endif

Kernels * Kernels_alloc(const char *fname);
double Kernels_Density(Kernels *kn, int ndim, double x[]);
void Kernels_free(Kernels * kn);

#ifdef __cplusplus
}
#endif

#endif
