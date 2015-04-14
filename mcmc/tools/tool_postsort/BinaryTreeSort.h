#ifndef BinaryTreeSort_h
#define BinaryTreeSort_h

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class Node
{
public:
	double val;
	int    index;
	Node * parent;
	Node * left;
	Node * right;

	//constructor
	Node()
	{
		parent = NULL;
		left   = NULL;
		right  = NULL;
	}
	Node(double _val, Node * _parent)
	{
		val = _val;
		parent = _parent;
		left   = NULL;
		right  = NULL;
	}
	Node(double _val, int _index, Node * _parent)
	{
		val = _val;
		index = _index;
		parent = _parent;
		left   = NULL;
		right  = NULL;
	}

	Node * firstchild();
};

class BinaryTreeSort
{
protected:
	Node * root;
	
	void Insert(double val);
	void Insert(double val, int index);

	// inorder trasversal
	Node * InorderTrasverse(Node * node01);

	// post order trasversal
	Node * PostorderTrasverse(Node * node01);

	// free tree
	void Free(void);
public:
	//constructor
	BinaryTreeSort()
	{
		root = NULL;
	}

	void Sort(vector<double>& array01, vector<double>& array02);
	void Sort(vector<double>& array, vector<int>& indice);
	void Sort(vector<double>& array);	
};

#endif
