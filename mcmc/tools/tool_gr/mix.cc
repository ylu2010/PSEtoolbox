#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <deque>
#include <set>
#include <cmath>

#include <values.h>

using namespace std;

vector<double> MixingFraction(deque< vector< vector<double> > >& data,
			      vector< unsigned char >& mask,
			      unsigned endpt, unsigned psize, bool verbose)
{
  if (verbose) {
				// Separator
    cout << setfill('-') << setw(72) << "-" << endl << setfill(' ');
  }
  


  //***************************************************
  //  Begin computation
  //***************************************************
  
				// Number of chains
  unsigned Mchain = data[0].size();
				// Number of values per chain 
				// (ignore mixture count)
  unsigned Nparam = data[0][0].size();


  unsigned begpt;
  if (endpt<psize) begpt = 0;
  else begpt = endpt-psize;
				// 
				// Count aberrant chains
				// 
  unsigned Mgood = 0;
  for (unsigned n=0; n<mask.size(); n++) 
    if (mask[n]==0) Mgood++;


				// Analyze chains
  vector<double> frac(Mchain, 1.0);
  for (unsigned n=0; n<Mchain; n++) {
    for (unsigned m=0; m<Nparam; m++) {
      set<double> s;
      for (unsigned k=begpt; k<endpt; k++) s.insert(data[k][n][m]);
      frac[n] = min<double>(frac[n], (double)s.size()/(endpt-begpt));
    }
  }
  
  if (verbose) {
    cout << "Fractions: ";
    unsigned origp = cout.precision(2);
    cout.setf(ios::fixed);
    for (unsigned n=0; n<Mchain; n++) cout << setw(5) << frac[n];
    cout << endl;
    cout.precision(origp);
    cout.unsetf(ios::fixed);
				// Separator
    cout << setfill('-') << setw(72) << "-" << endl << setfill(' ');
  }

  return frac;
}
