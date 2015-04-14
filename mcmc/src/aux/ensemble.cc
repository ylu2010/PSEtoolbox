#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "ensemble.H"

using namespace std;

Ensemble::Ensemble(const char ifname[])
{
	string line;
	istringstream sline;

	ifstream ifile(ifname);
	if(!ifile.is_open())
	{
		cout << "Ensemble: error in openning file." << endl;
		abort();
	}

	// load the file header and count how many columns there are
	getline(ifile, line);
	sline.str(line);
	ndim = 0;
	while(1)
	{
	    string token;
	    sline >> token;
	    if(token.size() == 0)
	        break;
	    ndim++;
	}
	ndim -= 1;
	sline.clear();

	// load the data
	while(1)
        {
	    double temp;
	    Particle point(ndim);

	    getline(ifile, line);
	    if(line.size() == 0)
	        break;
	    sline.str(line);
	    sline >> temp;	// log pdf
	    pdf.push_back(temp);
	    for(unsigned i=0; i<ndim; i++)
	    {
	        sline >> temp;
	        point.coord[i] = temp;
	    }
	    sline.clear();
	    particle.push_back(point);
        }
	ifile.close();
	ndim = ndim;
	nparticle = particle.size();

	cout << "Data loaded:" << endl;
	cout << "  There are " << nparticle << " points" << endl;
	cout << "  in " << ndim << " dimensional space." << endl;
}

// c interface
extern "C" {
Ensemble * Ensemble_alloc(const char ifname[])
{
        Ensemble * temp = new Ensemble(ifname);
        return temp;
}

void Ensemble_free(Ensemble * e)
{
        delete e;
}

}
