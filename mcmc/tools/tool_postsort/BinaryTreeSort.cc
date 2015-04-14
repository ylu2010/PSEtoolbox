#include <iostream>
#include <cmath>
#include <vector>

#include "BinaryTreeSort.h"

Node * Node::firstchild()
{
	if(left==NULL && right==NULL)
		return NULL;
	else if(left != NULL)
		return left;
	else
		return right;
}

void BinaryTreeSort::Insert(double val)
{
	Node * p=root;

	if(root == NULL)
		root = new Node(val, NULL);
	else
	{
		while(true)
		{
			if(val >= p->val)
			{
				if(p->right != NULL)
					p = p->right;
				else
				{
					p->right = new Node(val, p);
					break;
				}	
			}
			else
			{
				if(p->left != NULL)
					p = p->left;
				else
				{
					p->left  = new Node(val, p);
					break;
				}
			}
		}
	}
}

void BinaryTreeSort::Insert(double val, int index)
{
	Node * p=root;

	if(root == NULL)
		root = new Node(val, index, NULL);
	else
	{
		while(true)
		{
			if(val >= p->val)
			{
				if(p->right != NULL)
					p = p->right;
				else
				{
					p->right = new Node(val, index, p);
					break;
				}	
			}
			else
			{
				if(p->left != NULL)
					p = p->left;
				else
				{
					p->left  = new Node(val, index, p);
					break;
				}
			}
		}
	}
}

Node * BinaryTreeSort::InorderTrasverse(Node * node01)
{
	Node * node02=node01;
	
	if(node02->right != NULL)
	{
		node02 = node02->right;
		while(node02->left != NULL)
		{
			node02 = node02->left;
		}
	}
	else
	{
		while(node02->parent != NULL && node02 == node02->parent->right)
		{
			node02 = node02->parent;
		}
		node02 = node02->parent;
	}

	return node02;
}

Node * BinaryTreeSort::PostorderTrasverse(Node * node01)
{
	Node * parent = node01->parent;
		
	if(parent == NULL)	// it is root!
	{
		return NULL;
	}
	else
	{
		if(node01 == parent->right)
		{
			return parent;
		}
		else if(node01 == parent->left)
		{
			if(parent->right == NULL)
			{
				return parent;
			}
			else	// it has a sibling
			{
				Node * p = parent->right;
				while(p->firstchild() != NULL)
					p = p->firstchild();
				return p;
			}
		}
	}
}


void BinaryTreeSort::Free(void)
{
	Node * p01=root;
	Node * p02=root;

	while(p02 != NULL)
	{
		p01 = p02;
		p02 = p01->firstchild();
	}

	while(p01 != NULL)
	{
		p02 = p01;
		p01 = PostorderTrasverse(p02);
		delete p02;
	}
	root = NULL;
}

void BinaryTreeSort::Sort(vector<double>& array01, vector<double>& array02)
{
	int n = array01.size();
	Node * p01, * p02;
	int k=0;

	//construct the binary tree
	for(int i=0; i<n; i++)
	{
		Insert(array01[i]);
	}

	//trasverse the tree
	//find the leftmost node 
	p01 = root;
	p02 = root;
	while(p02 != NULL)
	{
		p01 = p02;
		p02 = p01->left;
	}
	k = 0;
	while(p01 != NULL)
	{
		array02[k] = p01->val;
		p01 = InorderTrasverse(p01);
		k++;	
	}

	//free tree
	Free();
}

void BinaryTreeSort::Sort(vector<double>& array, vector<int>& indice)
{
	int n = array.size();
	Node * p01, * p02;
	int k=0;

	//construct the binary tree
	for(int i=0; i<n; i++)
	{
		Insert(array[i], indice[i]);
	}

	//trasverse the tree
	//find the leftmost node 
	p01 = root;
	p02 = root;
	while(p02 != NULL)
	{
		p01 = p02;
		p02 = p01->left;
	}
	k = 0;
	while(p01 != NULL)
	{
		indice[k] = p01->index;
		p01 = InorderTrasverse(p01);
		k++;	
	}

	//free tree
	Free();
}

void BinaryTreeSort::Sort(vector<double>& array)
{
	int n = array.size();
	Node * p01, * p02;
	int k=0;

	//construct the binary tree
	for(int i=0; i<n; i++)
	{
		Insert(array[i]);
	}

	//trasverse the tree
	//find the leftmost node
	p01 = root;
	p02 = root;
	while(p02 != NULL)
	{
		p01 = p02;
		p02 = p01->left;
	}
	k = 0;
	while(p01 != NULL)
	{
		array[k] = p01->val;
		p01 = InorderTrasverse(p01);
		k++;	
	}

	// free tree
	Free();
}
