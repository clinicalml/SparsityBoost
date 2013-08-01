#include "data.h"

// -----------------------------------
int data::read(istream &in) { 
	assert(nvar>0); // this should be known ahead of time
	for (int i=0;i<nvar;i++) in >> x[i]; 
	return 1; 
}

// -----------------------------------
data::data() { 
	prev = NULL;
	nvar = 0; 
	x = NULL;
}

// -----------------------------------
data::data(int nvar) { 
	prev = NULL;
	this->nvar = nvar; assert(nvar>0);
	x = new int[nvar];
	for (int i=0;i<nvar;i++) x[i] = 0;
}

// -----------------------------------
data::~data() { 
	if (nvar>0) delete x; 
	nvar = 0; x = NULL;
	prev = NULL;
}