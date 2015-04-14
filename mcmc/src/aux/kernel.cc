#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "kernel.H"

using namespace std;

double Kernel::gaussian(vector<double> x)
{
	try {
		if(x.size() != coord.size())
			throw "inconsistency!";
	}
	catch(char * str) {
		cout << "Exception raised: " << str << endl;
	}

	double w = 0.;
	for(unsigned int i=0; i<ndim; i++)
		w += -(x[i]-coord[i])*(x[i]-coord[i]);
	w /= 2.*size*size;
	w -= 0.5*ndim*log(2.*M_PI) + ndim*log(size);

	return w;
}


Kernels::Kernels(const char * fname)
{
	unsigned int n_col;
	string line;
	istringstream sline;
	ifstream fin;

	fin.open(fname);
	try {
		if(!fin.is_open())
			throw "error in openning file.";
	}
	catch (string e) {
		cout << "Exception raised: " << e << endl;
	}
	
	// detect the number of columns
	getline(fin, line);
	sline.str(line);
	n_col = 0;
	while(1)
	{
		string token;
		sline >> token;
		if(token.size() == 0)
			break;
		n_col++;
	}
	sline.clear();
	ndim = n_col - 2;

	// load data
	while(1)
	{
		double t;
		Kernel temp;
		temp.ndim = ndim;
		temp.coord.assign(ndim, 0);

		getline(fin, line);
		if(line.size() == 0)
			break;

		sline.str(line);
		sline >> temp.pdf;
		sline >> temp.size;
		for(unsigned i=0; i<ndim; i++)
			sline >> temp.coord[i];
		kernel.push_back(temp);		
		sline.clear();
	}
	fin.close();

	// remove kernels with 0 size
	for(vector<Kernel>::iterator it = kernel.begin(); it != kernel.end(); )
	{
		if(it->size < 0.01)    // this number is subject to change
			kernel.erase(it);
		else
			it++;
	}	
	
	num_points = kernel.size();

	cout << "There are " << num_points << " in " 
	     << ndim << " dimensional space." << endl;
}

double Kernels::Density(vector<double> x)
{
	double pdf = 0.;
	for(unsigned i=0; i<num_points; i++) {
//		cout << exp( kernel[i].gaussian(x) ) << endl;
		pdf += exp( kernel[i].gaussian(x) );
	}
	pdf /= (double)num_points;

	return pdf;
}

extern "C" {

Kernels * Kernels_alloc(const char fname[])
{
	Kernels * kn = new Kernels(fname);
	return kn;
}

double Kernels_Density(Kernels * kn, int ndim, double x[])
{
	vector<double> _x;
	for(unsigned i=0; i<ndim; i++)
		_x.push_back(x[i]);
	return kn->Density(_x);
}

void Kernels_free(Kernels *kn)
{
	delete kn;
}

}

