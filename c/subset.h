#ifndef _subset_H
#define _subset_H

#include "utils.h"

class subset { 
public:
	int nodes, n; 
	int *ind; 
	
	void set(); 
	void set(subset *S); 
	void set(subset *S, int k); // skip k'th index
	int next(); 
	
	subset(int nodes, int n); 
	~subset();
};

class subset_array {
	float *A; 
	unsigned long **step; 
	unsigned long maxindex; 
public:
	int nodes, n; 

	void set(float val); 
	void set(subset *S, float val); 
	void set(int *I, int ni, float val); 

	float value(subset *S); 
	float value(int *I, int ni, int k); // skip k'th index
	
	subset_array(int nodes,int n); 
	~subset_array();
};

#endif