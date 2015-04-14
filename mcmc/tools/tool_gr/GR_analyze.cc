// -*- C++ -*-

/* 
   Reads in a multilevel output file and analyzes convergence
   using Gelman and Rubin (1992) with outlier detection and
   Bayes evidence computation.
*/


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>

using namespace std;

bool ComputeOutlierMask(vector<unsigned char>& mask,
			vector< vector<double> >& prb,
			vector< vector<double> >& z,
			double alpha, int maxoutlier, 
			double poffset, bool verbose);

double GR(deque< vector< vector<double> > >& data,
	  vector< unsigned char >& mask,
	  unsigned endpt, unsigned psize, 
	  bool verbose, bool debug);

vector<double> MixingFraction(deque< vector< vector<double> > >& data,
			      vector< unsigned char >& mask,
			      unsigned endpt, unsigned psize, bool verbose);


string inputfile("multilevel.dat");
string outputfile("new.chain");

/*
  Parse one line of the state log file

  Return true if level matches current level and populate data vectors

  Return false if level does not match and do not populate data vectors
*/
bool parse_line(char *line, int level, int& iter, vector<double>& prob,
		unsigned& m, vector<double>& parm, bool mixture)
{
  double v;
  int l;

  istringstream istr(line);	// check level
  istr >> l;
  if (l!=level) return false;

  istr >> iter;		// dump iteration counter

  prob.erase(prob.begin(), prob.end());

  istr >> v;			// read probability
  prob.push_back(v);
  istr >> v;			// read likeliehood
  prob.push_back(v);
  istr >> v;			// read prior
  prob.push_back(v);
  if (mixture) istr >> m;	// read number
  else m = 1;
				// read parameter vector
  parm.erase(parm.begin(), parm.end());
  while (istr) {
    istr >> v;
    if (!istr) break;
    parm.push_back(v);
  }

  return true;
}

//
// Parse the header into field labels.  Look for "Number" after
// the "Prior" field as a sign of a mixture
//
bool isMixture(char *line)
{
  vector<string> labs;
  char *p = line;
  while (*p!='\0') {
    if (*p=='"') {
      labs.push_back(string(""));
      while (*(++p)!='"') labs.back().push_back(*p);
    }
    ++p;
  }

  bool ret = false;
  for (unsigned n=0; n<labs.size(); n++) {
    if (labs[n].compare("Prior")==0) {
      if (labs[n+1].compare("Number")==0) ret = true;
      break;
    }
  }
  
  return ret;
}

int
main(int argc, char** argv)
{
  /***************************************************
    Global variables
   ***************************************************/
				// Use this level
  int level = 0;
				// Number of states in ensemble before testing
  int maxit = 500;
				// Number of steps between tests
  int nskip = 100;
				// Confidence interval for two-sided 
  double alpha = 0.05;		// Grubbs test

				// Maximum offset for accepting
  double poffset = -30.0;	// outlier in addition to Grubbs test


  int maxout = 6;		// Maximum number of outliers

  double MaxR = 1.2;		// GR statistic threshold

  bool dump = false;		// Dump converged, burnt-in chain

				// Retain last convergence test sample
  bool retain = true;		// in converged output

				// Compute Bayes Evidence using
  bool evidence = false;	// trimmed harmonic means


  /***************************************************
    Parse command line
   ***************************************************/

  int c;
  while (1) {
    c = getopt (argc, argv, "i:o:l:m:n:N:a:r:defh");
    if (c == -1) break;
     
    switch (c)
      {
      case 'l': level = atoi(optarg); break;
      case 'N': maxit = atoi(optarg); break;
      case 'n': nskip = atoi(optarg); break;
      case 'm': maxout = atoi(optarg); break;
      case 'a': alpha = atof(optarg); break;
      case 'r': MaxR = atof(optarg); break;
      case 'd': dump = true; break;
      case 'f': retain = false; break;
      case 'e': evidence = true; break;
      case 'i': inputfile.erase(); inputfile = optarg; break;
      case 'o': outputfile.erase(); outputfile = optarg; break;
      case 'h':
	cerr << endl << "**** GRanalyze *****" << endl << endl <<
	  "\tComputes the Gelman-Rubin statistic for DifferentialEvolution\n"
	  "\toutput.  Outlier chains are dropped automatically using the\n"
	  "\tGrubbs' statistic with the provided confidence region.  Mixing\n"
	  "\tfractions are also computed.  Poor mixing within chains will\n"
	  "\tcause Gelman-Rubin to always fail.\n\n";
      case '?': 
	string msg = "usage: " + string(argv[0]) + " [options]\n\n\twhere options are:\n\n";
	msg += "\t-l int   \t\tlevel to analyze (default: 0)\n";
	msg += "\t-N int   \t\tnumber of states in each analysis (default: 500)\n";
	msg += "\t-n int   \t\tanalysis interval (default: 100)\n";
	msg += "\t-m int   \t\tmaximum number of allowed outliers (default: 6)\n";
	msg += "\t-a float \t\tGrubbs' test confidence (default: 0.05)\n";
	msg += "\t-r float \t\tupper limit for GR statistic (default: 1.2)\n";
	msg += "\t-i string\t\tinput file (default: " + inputfile + ")\n";
	msg += "\t-d       \t\tdump converged chain with outliers removed\n";
	msg += "\t-o string\t\toutput file (default: " + outputfile + ")\n";
	msg += "\t-e       \t\tCompute Bayes evidence\n";
	msg += "\t-f       \t\tDO NOT dump tested steps with converged chain\n";
	cerr << msg; exit(-1);
      }
  }

  /***************************************************
    Read in multilevel data
   **************************************************/

  ifstream in(inputfile.c_str());
  if (!in) {
    cerr << argv[0] << ": error opening <" << inputfile << ">\n";
    exit(-1);
  }

  const int linesize=2048;
  char line[linesize], header[linesize];
  int iter, last=-1;
  unsigned m;
  vector<double> prb, parm;

				// Read in (and save) header string
  in.getline(header, linesize);
				// Sanity check
  if (string(header).find("Probability") == string::npos) {
    cerr << "Expect the first line to be a field header, but found:" << endl
	 << header << endl;
    exit(-1);
  }
  bool mixture = isMixture(header);

  deque< vector< vector<double> > > prob, data;
  deque< vector<unsigned> > numb;

				// Read in the lines at the desired level
  bool first = true;

  vector<unsigned> mvec;
  vector< vector<double> > prec, drec;

  while (in) {
    in.getline(line, linesize);
    if (!in) break;
    while (parse_line(line, level, iter, prb, m, parm, mixture)) {
      if (first || last!=iter) {
	if (mvec.size()) {
	  prob.push_back(prec);
	  numb.push_back(mvec);
	  data.push_back(drec);

	  prec.erase(prec.begin(), prec.end());
	  mvec.erase(mvec.begin(), mvec.end());
	  drec.erase(drec.begin(), drec.end());
	}
	last = iter;
	first = false;
      }
	
      prec.push_back(prb);
      mvec.push_back(m);
      drec.push_back(parm);

      in.getline(line, linesize);
      if (!in) break;
    }
    if (!first) break;
  }

  unsigned psize, dsize, rsize;

  if (prob.size() < 2 || data.size() < 2) {
    cout << "No states at Level " << level << endl;
    exit(-1);
  }

  cout << "=========================" << endl
       << "====  Data file info ====" << endl
       << "=========================" << endl
       << "Size of prob:       " << prob.size() << endl
       << "Size of data:       " << (dsize=data.size()) << endl
       << "Size of chain set:  " << (psize=data.front().size()) << endl
       << "Size of param vec:  " << (rsize=data.front()[0].size()) << endl;
  
  for (unsigned n=1; n<data.size(); n++) {
    if (data[n].size() != psize) {
      cout << "Record " << n << " has size=" << data[n].size() << endl;
    }
  }

  int nconverge=-1;
  vector<unsigned char> mask;
  unsigned nend=max<unsigned>(maxit, nskip), count=0, ssize=maxit;

  while (nend<dsize) {

    count++;

    if (maxit==0) ssize = nend/2;

    ostringstream ostr;
    ostr << "==== Interval #" << count << ": [" 
	 << max<int>(0, (int)nend-(int)maxit) << ", " << nend << "] ";
    cout << setw(72) << setfill('=') << "=" << endl
	 << setw(72) << left<< ostr.str() << endl
	 << setw(72) << setfill('=') << "=" << endl 
	 << setfill(' ');


    if (ComputeOutlierMask(mask, prob[nend], data[nend], 
			   alpha, maxout, poffset, true) )
      {
	vector<double> frac = MixingFraction(data, mask, nend, ssize, false);
	double wmix=2.0;
	unsigned wnum=-1;
	for (unsigned n=0; n<frac.size(); n++) {
	  if (mask[n]==1) continue;
	  if (frac[n]<wmix) {
	    wmix = frac[n];
	    wnum = n;
	  }
	}

	cout << "Worst non-masked mixing fraction in Chain #" << wnum
	     << " frac=" << wmix << endl;
	
	if (GR(data, mask, nend, ssize, true, false) < MaxR) {
	  if (retain)
	    nconverge = max<int>(0, nend-ssize);
	  else
	    nconverge = nend;
	  break;
	}
      }

    nend += nskip;
  }

  cout << endl;
  if (nconverge>=0)
    cout << "Converged at n=" << nconverge << endl;
  else
    cout << "No convergence" << endl;
    

  if (dump && nconverge>=0) {
    unsigned nchain = mask.size();
    unsigned nmask = 0;
    vector<unsigned char>::iterator ip;
    for (ip=mask.begin(); ip!=mask.end(); ip++) {
      if (*ip == 1) nmask++;
    }
    ofstream out(outputfile.c_str());
    out.precision(12);
    if (out) {
      cout << "Dumping converged states to file <"
	   << outputfile << "> with " << nmask << " chains removed" 
	   << endl;
      
      out << header << endl;
      for (unsigned i=nconverge; i<dsize; i++) {
	for (unsigned j=0; j<nchain; j++) {
	  if (mask[j]) continue;
	  out << setw(6) << level
	      << setw(7) << i;
	  for (int k=0; k<3; k++)
	    out << setw(20) << prob[i][j][k];
	  if (mixture) out << setw(10) << numb[i][j];
	  for (unsigned k=0; k<rsize; k++)
	    out << setw(20) << data[i][j][k];
	  out << endl;
	}
      }
    } else {
      cout << "Could not open dump file <" << outputfile << ">" << endl;
    }
  }

  if (evidence && nconverge>=0) {
    unsigned nchain = mask.size();
    unsigned nmask = 0;
    vector<unsigned char>::iterator ip;
    for (ip=mask.begin(); ip!=mask.end(); ip++) {
      if (*ip == 1) nmask++;
    }
    string file = outputfile + ".evidence";
    ofstream out(file.c_str());
    if (out) {
      cout << "Computing accumulated evidence with " << nmask 
	   << " chains removed" 
	   << endl;
      
      vector<double> llist;
      for (unsigned i=nconverge; i<dsize; i++) {
	for (unsigned j=0; j<nchain; j++) {
	  if (mask[j]) continue;
	  
	  llist.push_back(prob[i][j][1]);
	}
      }

      cout << "Using " << llist.size() << " states" << endl;
      sort(llist.begin(), llist.end(), greater<double>());

      double llmax = llist[0];
      double Lmax = exp(llist[0]);
      double kappa=1.0, low=0.0;

      for (unsigned t=1; t<llist.size(); t++) {
	kappa += exp(llmax - llist[t]);
	low += exp(llist[t] - llist[t-1]);
	out << setw(15) << llist[t]
	    << setw(15) << 1.0 - (1.0*t-1.0)/(llist.size()-1)
	    << setw(15) << kappa/t
	    << setw(15) << Lmax*low/kappa
	    << setw(15) << Lmax*t/kappa
	    << setw(15) << 0.5*Lmax*(low + t)/kappa
	    << endl;
      }
      
      cout << "Evidence=" 
	   << log(0.5*Lmax*(low + llist.size()-1)/kappa) << endl;

    } else {
      cout << "Could not open evidence file <" << file << ">" << endl;
    }
  }

  return 0;
}

