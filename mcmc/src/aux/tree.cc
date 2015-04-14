#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "ensemble.H"
#include "tree.H"

using namespace std;

unsigned KDTree::Insert(unsigned int index)
{
	bool same = false;

	// the tree hasnt been constructed
	if(root == NULL)
	{
		root = new Node();
		root->axis = 0;
		root->index = 0;
		root->coord = ensemble->particle[0].coord;
	}
	else // the tree is already there
	{
		Node * p = root;
		while(1)
		{
			// if it is the same point, don't do anything
			same = true;
			for(unsigned i=0; i<ndim; i++)
				same = same && (ensemble->particle[index].coord[i] == p->coord[i]);
			if(same)
				break;

			if(ensemble->particle[index].coord[p->axis] >= p->coord[p->axis])
			{
				if(p->childa == NULL)
				{
					p->childa = new Node();
					p->childa->axis = (p->axis+1)%ndim;
					p->childa->index = index;
					p->childa->coord = ensemble->particle[index].coord;
					p->childa->parent = p;
					break;
				}
				else
					p = p->childa;
			}
			else
			{
				if(p->childb == NULL)
				{
					p->childb = new Node();
					p->childb->axis = (p->axis+1)%ndim;
					p->childb->index = index;
					p->childb->coord = ensemble->particle[index].coord;
					p->childb->parent = p;
					break;
				}
				else
					p = p->childb;
			}
		}
	}

	if(same) return 0;
	else     return 1;
}

void KDTree::Construct()
{
	for(unsigned i=0; i<nparticle; i++)
	{
//		cout << "# " << i << endl;
		unsigned t;
		t = Insert(i);
		nnode += t;
		if(t)
			mask[i] = true;
		else
			mask[i] = false;
	}

	cout << "There are " << nnode << " nodes on the tree" << endl;
	cout << "with " << nparticle-nnode << " points removed." << endl;
}

int KDTree::Neighbor_WithinR(vector<double> coord, double r, vector<unsigned int> &indice)
{
	SearchTree(root, coord, r, indice);

	return 0;
}

int KDTree::Neighbor_kNearest(vector<double> coord, const unsigned int k, vector<double> &d, vector<unsigned int> &indice)
{
	try {	
		if(k >= nparticle)
			throw "Query not allowed.";
	}
	catch(char * str) {
		cout << "Exception raised: " << str << endl;
	}

	// init
	d.assign(k, 0.);
	indice.assign(k, 0);

	for(unsigned int i=0; i<k; i++)
	{
		indice[i] = i;
		d[i] = 0.;
		for(unsigned int j=0; j<ndim; j++)
			d[i] += (coord[j]-ensemble->particle[i].coord[j])*
			        (coord[j]-ensemble->particle[i].coord[j]);
		d[i] = sqrt(d[i]);
	}
	// put the maximum at the end
	unsigned int t=0;
	for(unsigned int it=0; it<k; it++)
	{
		if(d[it] > d[t])
			t = it;
	}
	swap(d[t], d[k-1]);
	swap(indice[t], indice[k-1]);

	// search the tree
	SearchTree(root, coord, d, indice);

	return 0;
}

int KDTree::Neighbor_Nearest(vector<double> coord, unsigned int &index)
{
	double d=0.0;

	// find the block this particle is in
	if(root == NULL)
		return -1;
	else
	{
		Node * p = root;
		while(1)
		{
			if(coord[p->axis] >= p->coord[p->axis])
			{
				if(p->childa == NULL)
				{
					index = p->index;
					for(int i=0; i<ndim; i++)
					{
						d += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
					}
					d = sqrt(d);
					break;
				}
				else
					p = p->childa;
			}
			else
			{
				if(p->childb == NULL)
				{
					index = p->index;
					for(int i=0; i<ndim; i++)
					{
						d += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
					}
					d = sqrt(d);
					break;
				}
				else
					p = p->childb;
			}
		}
	}

	SearchTree(root, coord, d, index);

	return 0;
}


int KDTree::Neighbor_Nearest(vector<double> coord, const unsigned id, unsigned int &index)
{
	double d=0.0;

	// pick up a non-id particle
	Node * p = root;
	while( p!= NULL)
	{
		if(p->index != id)
			break;
		p = Preorder(p);
	}
	for(unsigned i=0; i<ndim; i++)
		d += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
	d = sqrt(d);

	SearchTree(root, coord, id, d, index);

	return 0;
}


int KDTree::Neighbor_Nearest(vector<Pair> &pairs)
{
	Node * p = root;
	Node * neighbor;
	double d;
	unsigned int index;

	// there must be at least 2 particles
	if(p == NULL || p->firstchild() == NULL)
		return -1;
	// walk the tree
	while(p != NULL)
	{
		// initialize d.
		if(p == root)
			neighbor = p->firstchild();
		else
			neighbor = p->parent;
		d = 0.0;
		for(int i=0; i<ndim; i++)
			d += (p->coord[i] - neighbor->coord[i]) * (p->coord[i] - neighbor->coord[i]);
		d = sqrt(d);
		index = neighbor->index;
		// search for the nearest neighbor
		SearchTree(root, p->coord, p->index, d, index);
		pairs[p->index].a = p->index;
		pairs[p->index].b = index;
		// move on
		p = Preorder(p);
	}

	return 0;
}

int KDTree::SearchTree(Node *p, const vector<double> coord, const double d, vector<unsigned int> &indice)
{
//	static unsigned int count = 0;
//	count++;
//	cout << "count=" << count << endl;

	if(p == NULL)
		return 0;
	else
	{
		double normal = coord[p->axis] - p->coord[p->axis];
		// first child
		if(normal >= 0.0 || fabs(normal) <= d)
			SearchTree(p->childa, coord, d, indice);
		// second child
		if(normal <= 0.0 || fabs(normal) <= d)
			SearchTree(p->childb, coord, d, indice);
		// parent
		double d_t = 0.0;
		for(int i=0; i<ndim; i++)
		{
			d_t += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
		}
		d_t = sqrt(d_t);
		if(d_t < d)
			indice.push_back(p->index);
		return 0;
	}
}

int KDTree::SearchTree(Node *p, const vector<double> coord, double &d, unsigned int &index)
//	recursive scheme
{
//	static unsigned int count = 0;
//	count++;
//	cout << "count=" << count << endl;

	if(p == NULL)
		return 0;
	else
	{
		double normal = coord[p->axis] - p->coord[p->axis];
		// first child
		if(normal >= 0.0 || fabs(normal) <= d)
			SearchTree(p->childa, coord, d, index);
		// second child
		if(normal <= 0.0 || fabs(normal) <= d)
			SearchTree(p->childb, coord, d, index);
		// parent
		double d_t = 0.0;
		for(int i=0; i<ndim; i++)
		{
			d_t += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
		}
		d_t = sqrt(d_t);
		if(d_t < d)
		{
			d = d_t;
			index = p->index;
		}
		return 0;
	}	
}

int KDTree::SearchTree(Node *p, const vector<double> coord, vector<double> &d, vector<unsigned int> &index)
//	recursive scheme
//	to find the k nearest points
{
//	static unsigned int count = 0;
//	count++;
//	cout << "count=" << count << endl;

	if(p == NULL)
		return 0;
	else
	{
		unsigned int size = d.size();
		double normal = coord[p->axis] - p->coord[p->axis];
		// first child
		if(normal >= 0.0 || fabs(normal) <= d[size-1])
			SearchTree(p->childa, coord, d, index);
		// second child
		if(normal <= 0.0 || fabs(normal) <= d[size-1])
			SearchTree(p->childb, coord, d, index);
		// parent
		double d_t = 0.0;
		for(int i=0; i<ndim; i++)
		{
			d_t += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
		}
		d_t = sqrt(d_t);
		if(d_t < d[size-1])
		{
			d[size-1] = d_t;
			index[size-1] = p->index;
		}

		// put the maximum at the end of vector d
		unsigned int t = 0;
		for(unsigned int it=0; it<size; it++)
		{
			if(d[it] > d[t])
				t = it;
		}
		swap(d[t], d[size-1]);
		swap(index[t], index[size-1]);

		return 0;
	}	
}

int KDTree::SearchTree(Node *p, const vector<double> coord, const unsigned int id, double &d, unsigned int &index)
//	recursive scheme
{
//	static unsigned int count = 0;
//	count++;
//	cout << "count=" << count << endl;

	if(p == NULL)
		return 0;
	else
	{
		double normal = coord[p->axis] - p->coord[p->axis];
		// first child
		if(normal >= 0.0 || fabs(normal) <= d)
			SearchTree(p->childa, coord, id, d, index);
		// second child
		if(normal <= 0.0 || fabs(normal) <= d)
			SearchTree(p->childb, coord, id, d, index);
		// parent
		if(p->index != id)
		{
			double d_t = 0.0;
			for(int i=0; i<ndim; i++)
			{
				d_t += (coord[i]-p->coord[i])*(coord[i]-p->coord[i]);
			}
			d_t = sqrt(d_t);
			if(d_t < d)
			{
				d = d_t;
				index = p->index;
			}
		}
		return 0;
	}	
}

Node * KDTree::Preorder(Node *p)
{
	if(p->firstchild() != NULL)
		return p->firstchild();
	else		// no child at all
	{
		Node * next = p;
		while(next != NULL)
		{
			if(next->sibling() == NULL)
				next = next->parent;
			else
				break;
		}
		if( next == NULL)
			return NULL;
		if( next->sibling() != NULL)
			return next->sibling();
	}
}

void KDTree::FreeTree()
{
	FreeNode(root);
}

void KDTree::FreeNode(Node *p)
{
	if(p != NULL)
	{
		FreeNode(p->childa);
		FreeNode(p->childb);
		delete p;
	}
}

// c interface
extern "C" {

KDTree * KDTree_alloc(Ensemble * e)
{
        KDTree * temp = new KDTree(e);
        return temp;
}

void KDTree_construct(KDTree *kd)
{
        kd->Construct();
}

double KDTree_getPDF(KDTree *kd, int ndim, double coord[])
{
	vector<double> _coord;
	for(unsigned i=0; i<ndim; i++)
		_coord.push_back(coord[i]);
        return kd->getPDF(_coord);
}

void KDTree_free(KDTree *kd)
{
        delete kd;
}

}
