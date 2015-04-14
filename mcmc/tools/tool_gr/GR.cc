#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <deque>
#include <cmath>

#include <values.h>

using namespace std;

double GR(deque< vector< vector<double> > >& data,
	  vector< unsigned char >& mask,
	  unsigned endpt, unsigned psize, 
	  bool verbose, bool debug)
{

  if (verbose) {
				// Separator
    cout << setfill('-') << setw(72) << "-" << endl
	 << "End point=" << endpt << "  Count=" << psize << endl
	 << setw(72) << "-" << endl << setfill(' ');
  }
  


  //***************************************************
  //  Begin computation
  //***************************************************
  
				// Number of chains
  unsigned Mchain = data[0].size();
				// Number of values per chain 
				// (ignore mixture count)
  unsigned Nparam = data[0][0].size();

				// Number of steps to use
  unsigned Nsteps = psize;

  unsigned begpt;
  if (endpt+1<psize) begpt = 0;
  else begpt = endpt+1-psize;
				// 
				// Count aberrant chains
				// 
  unsigned Mgood = 0;
  for (unsigned n=0; n<mask.size(); n++) 
    if (mask[n]==0) Mgood++;

				// 
				// Gelman & Rubin quantities
				// 
  vector<double> B(Nparam, 0.0), W(Nparam, 0.0), meanS(Nparam, 0.0);
  vector<double> Sig2Hat(Nparam, 0.0), rHat(Nparam, 0.0);

				// 
				// Set up vectors for in-chain analysis
				// 
  vector< vector<double> > meanCh(Mchain), varCh(Mchain);
  for (unsigned i=0; i<Mchain; i++) {
    if (mask[i]) continue;
    meanCh[i]  = vector<double>(Nparam, 0.0);
    varCh[i]   = vector<double>(Nparam, 0.0);
  }
				// 
				// Compute means and variance of each chain
				//
  for (unsigned n=begpt; n<=endpt; n++) {
    for (unsigned i=0; i<Mchain; i++) {
      if (mask[i]) continue;
      for (unsigned j=0; j<Nparam; j++) {
	meanCh[i][j]  += data[n][i][j];
	varCh [i][j]  += data[n][i][j]*data[n][i][j];
      }
    }
  }

  for (unsigned i=0; i<Mchain; i++) {
    if (mask[i]) continue;
    for (unsigned j=0; j<Nparam; j++) {
      meanCh[i][j] /= Nsteps;
      varCh[i][j] = 
	(varCh[i][j] - meanCh[i][j]*meanCh[i][j]*Nsteps)/(Nsteps-1);
      meanS[j] += meanCh[i][j]/Mgood;
    }
  }
				// 
				// Compute the Gelman & Rubin quantities
				// 
  
  for (unsigned i=0; i<Mchain; i++) {
    if (mask[i]) continue;
    for (unsigned j=0; j<Nparam; j++) {
      B[j] += 
	(meanCh[i][j] - meanS[j])*(meanCh[i][j] - meanS[j]) *
	Nsteps/(Mgood-1);
      W[j] += varCh[i][j]/Mgood;
    }
  }
    
  
  for (unsigned j=0; j<Nparam; j++) {
    if (W[j]>0.0) {
      Sig2Hat[j] = (W[j]*(Nsteps-1) + B[j])/Nsteps;
      rHat[j] = sqrt(Sig2Hat[j]/W[j]);
    }
  }
  
  
  if (verbose && debug) {
    cout << endl 
	 << "Per chain data, N=" << Nsteps
	 << " M=" << Mgood << "/" << Mchain
	 << ": " << endl << endl
	 << setw(4) << left << "#";
    for (unsigned j=0; j<Nparam; j++) {
      ostringstream out1, out2;
      out1 << "Mean [" << j << "]";
      out2 << "Var [" << j << "]";
      cout << setw(15) << left << out1.str()
	   << setw(15) << left << out2.str();
    }
    cout << endl
	 << setw(4) << left << "-";
    for (unsigned j=0; j<Nparam; j++) {
      ostringstream out1, out2;
      out1 << "Mean [" << j << "]";
      out2 << "Var [" << j << "]";
      cout << setw(15) << left << "---------"
	   << setw(15) << left << "---------";
    }
    cout << endl;

    for (unsigned i=0; i<Mchain; i++) {
      cout << setw(4) << i;
      if (mask[i])
	for (unsigned j=0; j<Nparam; j++) {
	  cout << setw(15) << left << "*******"
	       << setw(15) << left << "*******";
	}
      else
	for (unsigned j=0; j<Nparam; j++) {
	  cout << setw(15) << left << meanCh[i][j]
	       << setw(15) << left << varCh[i][j];
	}
      cout << endl;
    }
    cout << endl << setw(4) << "**";
    for (unsigned j=0; j<Nparam; j++) {
      cout << setw(15) << left << meanS[j]
	   << setw(15) << left << B[j]/Nsteps;
    }
    cout << endl;
  }
      

  if (verbose) cout << endl 
		    << "Convergence summary, N=" << Nsteps
		    << " M=" << Mgood << "/" << Mchain
		    << ": " << endl << endl
		    << setw(4) << left << "#"
		    << setw(15) << left << "Mean"
		    << setw(15) << left << "Between"
		    << setw(15) << left << "Within"
		    << setw(15) << left << "Total"
		    << setw(15) << left << "Rhat"
		    << endl
		    << setw(4) << left << "-"
		    << setw(15) << left << "----"
		    << setw(15) << left << "-------"
		    << setw(15) << left << "------"
		    << setw(15) << left << "-----"
		    << setw(15) << left << "----"
		    << endl;
  
  double minParam = MAXDOUBLE, maxParam = 0.0;
  for (unsigned j=0; j<Nparam; j++) {

    if (verbose) cout << setw(4)  << j
		      << setw(15) << meanS[j]
		      << setw(15) << B[j]/Nsteps
		      << setw(15) << W[j]
		      << setw(15) << Sig2Hat[j]
		      << setw(15) << rHat[j]
		      << endl;

    minParam = min<double>(rHat[j], minParam);
    maxParam = max<double>(rHat[j], maxParam);
  }

  if (verbose) {
    cout << endl
	 << "Min, Max = " << minParam << ", " << maxParam << endl;
  }

  return maxParam;
}
