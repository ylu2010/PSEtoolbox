#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>

#include "state.h"
#include "BinaryTreeSort.h"

using namespace std;

int main(int argc, char * argv[])
{
	string ifname("turtle.dat");
	string ofname("turtle.dat.top");
	double p=1.0;
	int ndim = 10;

	while(1)
	{
		int c = getopt(argc, argv, "i:o:p:d:");
		if( c == -1 ) break;
		switch(c)
		{
			case 'i': ifname   = optarg; 		break;
			case 'o': ofname   = optarg;		break;
			case 'p': p	   = atof(optarg);	break;
			case 'd': ndim	   = atoi(optarg);	break;
			case '?': 
				cerr << "please offer a value for that each option." << endl;
				abort();
		}	
	}
	
	vector<double> posteriers;
	double	       posterier;
	string         temp;
	vector<int>    indice;
	int            index = 0;
	vector<string> lines;
	string         line;
	istringstream  lline;
	string         header;

	BinaryTreeSort sort;

	int nstate;
	int cutoff;

	ifstream ifile(ifname.c_str());
	ofstream ofile(ofname.c_str());

//	load original file
	getline(ifile, header);
	while(!ifile.eof())
	{
		getline(ifile, line);
		lines.push_back(line);
		lline.str(line);		
		lline >> temp;		// level
		lline >> temp;		// iter
		lline >> posterier;	// posterier
		posteriers.push_back(posterier);
		indice.push_back(index);
		index++;
		lline.clear();
	}
	ifile.close();

//	sort
	sort.Sort(posteriers, indice);

//	output
	nstate = indice.size();
	cutoff = nstate*(1.-p);
	ofile << header;
	for(int i=cutoff; i<nstate; i++)
	{
		ofile << endl;
		ofile << lines[indice[i]];
	}
	ofile.close();

//	end
	return 0;
}
