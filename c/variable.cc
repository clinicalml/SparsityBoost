#include "variable.h"

// ---------------------------------------------
parentset::parentset() { 
	npa = 0; pa = NULL; 
	W = 0; prev = NULL;
}

// ---------------------------------------------
parentset::parentset(int npa) { 
	this->npa = npa;
	if (npa>0) pa = new int[npa]; else pa = NULL; 
	W = 0; 	prev = NULL;
}

// ---------------------------------------------
parentset::parentset(int npa, int *pa, float W) { 
	this->npa = npa;
	if (npa>0) this->pa = new int[npa]; else this->pa = NULL; 
	for (int j=0;j<npa;j++) this->pa[j] = pa[j]; 
	this->W = W; prev = NULL;
}

// ---------------------------------------------
parentset::~parentset() { 
	if (npa>0) delete pa;
	npa = 0; pa = NULL; 
	W = 0; 	
	prev = NULL;
}

// ---------------------------------------------
void delete_parentsets(parentset *paset) { 
	for (parentset *s = paset; s != NULL; ) {
		parentset *tmp = s; s = s->prev; delete tmp;
	}
}

// ---------------------------------------------
void parentset::read(istream &in) { 
	in >> W >> npa; assert(npa>=0); 
	if (npa>0) { 
		pa = new int[npa]; 
		for (int j=0;j<npa;j++) in >> pa[j]; 
	} else pa = NULL;
}

// ---------------------------------------------
void parentset::print(ostream &out) { 
	out << W << " " << npa << " ";
	for (int j=0;j<npa;j++) out << pa[j] << " ";
	out << endl; 
}

// =============================================================================

void variable::print(ostream &out) {
	unsigned long npasets = 0; 
	for (parentset *s = paset; s != NULL; s = s->prev) npasets++;
	
	out << id << " " << npasets << endl;
	for (parentset *s = paset; s != NULL; s = s->prev) s->print(out);
}

// ---------------------------------------------
variable::variable() {
	id = 0; paset = NULL;
	prev = NULL;
}

// ---------------------------------------------
void delete_variables(variable *var) { 
	for (variable *v = var; v != NULL; ) {
		variable *tmp = v; v = v->prev; delete tmp;
	}
}

// ---------------------------------------------
variable::~variable() {
	if (paset != NULL) { 
		delete_parentsets(paset); 
	}
	id = 0; paset = NULL;
	prev = NULL;
}
