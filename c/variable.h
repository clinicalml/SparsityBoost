#ifndef _variable_H
#define _variable_H

#include "utils.h"

class parentset;
class parentset { 
public:
	int npa, *pa; 
	float W; 
	parentset *prev; 

	void read(istream &in); 
	void print(ostream &out); 
	
	parentset(); 
	parentset(int npa); 
	parentset(int npa, int *pa, float W); 
	~parentset();
};
void delete_parentsets(parentset* paset); 

class variable;
class variable { 
public:	
	int id; 
	parentset *paset; 
	variable *prev; 

	void print(ostream& out); 
	
	variable();
	~variable();
};

void delete_variables(variable *var); 
#endif
