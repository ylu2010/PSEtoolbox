#ifndef ensemble_h
#define ensemble_h

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Particle
{
public:
	unsigned int ndim;
	vector<double> coord;

	// constructor
	Particle(unsigned int _ndim)
	{
		ndim = _ndim;
		coord.assign(_ndim, 0.0);
	}

	// Distance between two particles
	double operator-(Particle b)
	{
		double d = 0.0;

		for(unsigned i=0; i<ndim; i++)
			d += (coord[i] - b.coord[i])*(coord[i] - b.coord[i]);
		d = sqrt(d);

		return d;
	}
};

class Ensemble
{
public:
	// # of particles
	unsigned int nparticle;
	// # of dimensionality
	unsigned int ndim;

	// particles and pdf
	vector<Particle> particle;
	vector<double> pdf;

	// constructor
	Ensemble() {}
	Ensemble(unsigned int _nparticle, unsigned int _ndim)
	{
		nparticle = _nparticle;
		ndim      = _ndim;
	}
	Ensemble(const char ifname[]);
	// print
//	printf(const char ofname[]);
};

#else
typedef
  struct Ensemble
  Ensemble;
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern Ensemble * Ensemble_alloc(const char ifname[]);
extern void Ensemble_free(Ensemble * e);

#ifdef __cplusplus
}
#endif

#endif
