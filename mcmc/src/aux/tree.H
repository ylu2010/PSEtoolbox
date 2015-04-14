#ifndef tree_h
#define tree_h

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <cstdlib>

#include "ensemble.H"

using namespace std;

class Node
{
public:
	unsigned int axis;
	int index;
	vector<double> coord;

	Node * parent;
	Node * childa;
	Node * childb;

	// default constructor
	Node()
	{
		parent = NULL;
		childa = NULL;
		childb = NULL;
	}

	// sibling
	Node * sibling()
	{
		if(parent == NULL)	// this == root
			return NULL;
		else if(this == parent->childa)
			return parent->childb;
		else if(this == parent->childb)
			return NULL;
		else 			// something must be wrong with the tree
			abort();
	}

	// firstchild
	Node * firstchild()
	{
		if(childa != NULL)
			return childa;
		else if(childb != NULL)
			return childb;
		else
			return NULL;
	}
};

class Pair
{
public:
	unsigned int a;
	unsigned int b;
};

class KDTree
{
protected:
	Ensemble * ensemble;

	unsigned int ndim;
	unsigned int nparticle;
	unsigned int nnode;
	vector<bool> mask;
public:
	Node * root;

	// constructor
	KDTree(Ensemble * _ensemble)
	{
		root = NULL;
		ensemble = _ensemble;
		ndim = ensemble->ndim;
		nparticle = ensemble->nparticle;
		nnode = 0;
		mask.assign(nparticle, true);
	}
	// destructor
	~KDTree()
	{
		FreeTree();
	}

	// construct the tree
	unsigned Insert(unsigned int index);
	void Construct();

	// access
	bool Mask(unsigned i)
	{
		return mask[i];
	}

	// the k nearest neighbor
	int Neighbor_kNearest(vector<double> coord, const unsigned k, vector<double> &d, vector<unsigned int> &indice);
	// the nearest neighbor
	int Neighbor_Nearest(vector<double> coord, unsigned int &index);
	// the nearest neighbor but skip the point with index == id
	int Neighbor_Nearest(vector<double> coord, const unsigned id, unsigned int &index);
	// nearest neighbor of all the particles
	int Neighbor_Nearest(vector<Pair> &pairs);
	// neighbors within radius r
	int Neighbor_WithinR(vector<double> coord, double r, vector<unsigned int> &indice);

	// to find the k nearest neighbor
	int SearchTree(Node *p, const vector<double> coord, vector<double> &d, vector<unsigned int> &index);
	// to find the nearest neighbor
	int SearchTree(Node *p, const vector<double> coord, double &d,      unsigned int &index);
	// return all the id's of the particles within radius d
	int SearchTree(Node *p, const vector<double> coord, const double d, vector<unsigned int> &indice);
	// skip the particle with index == id
	int SearchTree(Node *p, const vector<double> coord, const unsigned int id, double &d, unsigned int &index);

	// trasverse the tree
	Node * Preorder(Node * p);

	// 
	double getPDF(vector<double> coord)
	{
		unsigned int index;
		Neighbor_Nearest(coord, index);
		return ensemble->pdf[index];
	}

	// destruct the tree
	void FreeTree();
	void FreeNode(Node *p);
};

#else

#include "ensemble.H"

typedef 
  struct KDTree
  KDTree;
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern KDTree * KDTree_alloc(Ensemble *e);
extern void KDTree_construct(KDTree *kd);
extern double KDTree_getPDF(KDTree *kd, int ndim, double coord[]);
extern void KDTree_free(KDTree *kd);

#ifdef __cplusplus
}
#endif

#endif
