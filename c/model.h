#ifndef _model_H
#define _model_H

#include "utils.h"
#include "variable.h"
#include "data.h"
#include "subset.h"
#include "LgBetaTable.h"
#include <limits>

typedef std::numeric_limits< double > dbl;

class model;
class model { 
	double score(subset*S,double alpha); 
	double bic_score(subset*S); 
	void edge_scores(subset* S, double** min_edge_scores, /*PyObject *pInstance,*/ LgBetaTable&);
	void edge_taus(subset* S, double** min_edge_scores);

public:
	int nodes; 
	int *nval;
	variable **VAR; 
	data *D; 
	double eta;
	
	void read_data(options& opt);
	void print_model(options& opt);
	void learn_model(options& opt);
	void calculate_epsilon(options& opt);
	
	model(); 
	~model();
};

#endif
