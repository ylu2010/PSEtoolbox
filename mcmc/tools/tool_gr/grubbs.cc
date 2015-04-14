#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>

using namespace std;

extern double inv_student_t_1sided(double alpha, double df);

/*
  Grubbs' for outlier detection.

  Grubbs' test (Grubbs 1969 and Stefansky 1972) is used to detect
  outliers in a univariate data set. It is based on the assumption of
  normality.  This should be the case for converged MCMC calc.
*/

bool ComputeOutlierMask(vector<unsigned char>& mask,
			vector< vector<double> >& prb,
			vector< vector<double> >& z,
			double alpha, int maxoutlier, 
			double poffset, bool verbose)
{


  unsigned n = z.size();
  unsigned m = z.front().size();
  int cnt=0;

  mask = vector<unsigned char>(n, 0);


  double maxLP = -1.0e20;
  for (unsigned i=0; i<n; i++) maxLP = max<double>(maxLP, prb[i][0]);

  //
  // Do each variable separately
  //
  for (unsigned q=0; q<m; q++) {
				// P is the sorted list of (value, index) pairs
    vector< pair<double, int> > p;
    for (unsigned i=0; i<n; i++) {
      if (mask[i]==0)
	p.push_back(pair<double, int>(z[i][q], i));
    }
  
				// Compute the Student statistic for 
				// k-cnt possible outliers
    double val, pcrit, t, G, gmax;
    unsigned sz=p.size();
    unsigned k = sz*3/4;
    if (k<2) break;

    double mean, var;

    for (unsigned j=0; j<k; j++) {
      mean = var = 0.0;		// Compute mean and variance 
				// for the remaining sample
      for (unsigned i=0; i<sz; i++) {
	val = z[p[i].second][q];
	mean += val/sz;
	var  += val*val/(sz-1);
      }
      var -= mean*mean*sz/(sz-1);

      if (var<1.0e-18) break;

      for (unsigned i=0; i<sz; i++) p[i].first = fabs(z[p[i].second][q]-mean);
      sort(p.begin(), p.end());

				// The supremum value
      G = p.back().first / sqrt(var);

				// 2-sided Grubb's test to find outliers
      pcrit = 1.0 - (alpha*0.5/sz);
    
      t = inv_student_t_1sided(pcrit, sz-2);

      gmax = t*(sz-1)/sqrt(((t*t + (sz-2))*sz));

				// Report
      if (verbose) {
	if (G > gmax && prb[p.back().second][0]+poffset<maxLP)
	  cout << "***** Outlier #" << cnt+1
	       << " in Chain #" << p.back().second 
	       << ", analysis:" << endl
	       << "  " << setw(20) << "Variable #"  << " = " << setw(15) << q << endl
	       << "  " << setw(20) << "Sample size" << " = " << setw(15) << sz << endl
	       << "  " << setw(20) << "Minimum"     << " = " << setw(15) << p[0].first << endl
	       << "  " << setw(20) << "Maximum"     << " = " << setw(15) << p[sz-1].first << endl
	       << "  " << setw(20) << "Mean"        << " = " << setw(15) << mean << endl
	       << "  " << setw(20) << "Std dev"     << " = " << setw(15) << sqrt(var) << endl
	       << "  " << setw(20) << "Grubbs test" << " = " << setw(15) << G << endl
	       << "  " << setw(20) << "Lambda"      << " = " << setw(15) << gmax << endl
	       << "  " << setw(20) << "Crit value"  << " = " << setw(15) << alpha << endl
	       << "  " << setw(20) << "logP offset" << " = " << setw(15) << prb[p.back().second][0]-maxLP << endl
	       << endl;
      }
      
				// Record the outlier
      if (G > gmax && prb[p.back().second][0]+poffset<maxLP) {
	mask[p.back().second] = 1;
	p.pop_back();		// and delete it from the list
	cnt++;
      }
      else break;			// No more outliers . . . we're done
    }
  }

  if (cnt > maxoutlier) {
				// Warning for the naive.
				// (If it's good enough for Braak 2006 it's
				//  good enough for me . . .)
    cout << "****** Too many aberrant chains!  ******" << endl
	<< "****** Will continue sans masking ******" << endl
	<< setw(72) << setfill('-') << "-" << endl << setfill(' ');

    return false;
  }

				// Set the new mask
  if (verbose) {
    if (cnt==0) cout << "No outliers" << endl;
    else {
      cout << "Outlier mask: ";
      for (unsigned i=0; i<n; i++) 
	if (mask[i]==0) cout << setw(2) << 0;
	else                cout << setw(2) << 1;
      cout << endl;
    }
    cout << setw(72) << setfill('-') << "-" << endl << setfill(' ');
  }

  return true;
}

