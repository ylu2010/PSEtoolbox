#ifndef state_h
#define state_h

#include <vector>
#include <iostream>

using namespace std;

class state{
public:
	int level;
	int iter;
	double posterier;
	double likelihood;
	double prior;

	int nparam;
	vector<double> param;

	// Constructor
	state(double _nparam)
	{
		nparam = _nparam;
		param = vector<double>(nparam, 0.0);
	}
	// Destructor
	// default	
};

#endif
