#include "model.h"
#include <cstring>
#include <algorithm>
//#include <Python.h>
#include <vector>
#define LEN 4096

#define MAX_CONDITIONING_SET_SIZE 2

// This is the weighting value, i.e. the amount that the data-driven complexity penalty should be weighted (\psi_2 in Eliot's writeup).
#define PSI2 1

// ---------------------------------------------------
double model::score(subset* S, double alpha) { 

	// prepare to index configurations
	int step[S->n+1]; int maxindx = 1; 
	for (int i=0;i<S->n;i++) { step[i]=maxindx; maxindx *= nval[S->ind[i]]; }
	
	double N[maxindx]; for (int i=0;i<maxindx;i++) N[i] = alpha/(double)maxindx; 

	// calculate the score
	double val = 0; 
	for (data *d=D;d != NULL;d=d->prev) { 
		int xindx = 0; for (int i=0;i<S->n;i++) xindx += d->x[S->ind[i]]*step[i]; 
		val += log(N[xindx]); N[xindx] += 1.0;
	}
	
	return val; 
}

// ---------------------------------------------------
// This function computes H(X_s) - .5 log(ndata)#assignments(S)
double model::bic_score(subset* S) { 
	
	// prepare to index configurations
	int step[S->n+1]; int maxindx = 1; 
	for (int i=0;i<S->n;i++) { step[i]=maxindx; maxindx *= nval[S->ind[i]]; }
	
	double N[maxindx]; for (int i=0;i<maxindx;i++) N[i] = 0; 
	
	// calculate the score
	int ndata = 0;
	for (data *d=D;d != NULL;d=d->prev) { 
		int xindx = 0; for (int i=0;i<S->n;i++) xindx += d->x[S->ind[i]]*step[i]; 
		N[xindx] += 1.0; ndata++; 
	}
	double val = 0; 
	for (int i=0;i<maxindx;i++) if (N[i]>0) val += N[i]*log(N[i]); 
	
	val += -0.5*maxindx*log(ndata); 
	
	return val; 
}


// ---------------------------------------------------
void model::edge_scores(subset* S, double** min_edge_scores, /* PyObject *pInstance,*/ LgBetaTable &etasLgBetaTable ) { 
	
  //PyObject *pValue;

  //  int total_beta_calls=0;
  //  int stopped_early_beta_calls=0;
  
  vector<std::pair<double, int> > MI_N_pairs;
  int total; double mi; // Used repeatedly below.
  double log_of_half = log(.5);

	// prepare to index configurations
	int step[S->n+1]; int maxindx = 1; 
	for (int i=0;i<S->n;i++) { step[i]=maxindx; maxindx *= nval[S->ind[i]]; }
	
	double N[maxindx]; for (int i=0;i<maxindx;i++) N[i] = 0; 
	//	double Ncopy[maxindx];
	
	// go through all of the data, computing counts for this subset
	int ndata = 0;
	for (data *d=D;d != NULL;d=d->prev) { 
		int xindx = 0; for (int i=0;i<S->n;i++) xindx += d->x[S->ind[i]]*step[i]; 
		N[xindx] += 1.0; ndata++; 
	}
	//	cout << "Number of lines:"  << ndata << endl;

	// for all pairs i & j and for each setting of the other variables, compute beta.
	int assignment[S->n];
	for (int vari=0; vari < S->n; vari++) {
	  for (int varj= vari+1; varj < S->n; varj++) {

	    /* Do this in two steps:
	       First compute all MI,N pairs. Then, sort by MI and evaluate beta. The reason we do this
	       is so that we can increase the probability that we can quit early due to pruning.
	       There will be k^subsetsize: can precompute if necessary (so no memory allocation)
	    */
	    MI_N_pairs.clear();

	    // Output this set
	    //	    out << S->ind[vari] << "," << S->ind[varj] << " ";
	    //	    for (int i=0; i<=S->n-1; i++)
	    //	      if(i!= vari && i != varj)
	    //		out << S->ind[i] << ",";
	    //	    out << " ";

	    double pairwise[nval[S->ind[vari]]][nval[S->ind[varj]]];
	    double marginal_vari[nval[S->ind[vari]]];
	    double marginal_varj[nval[S->ind[varj]]];

	    // initialize assignment to all zeros
	    for (int i=0; i<S->n; i++) assignment[i] = 0;

	    // for each value of the conditioning set
	    bool done = false;
	    while(!done) {

	      // compute the pairwise distribution
	      int indx_start = 0;
	      for (int i=0;i<S->n;i++)
		if(i!= vari && i != varj)
		  indx_start += assignment[i]*step[i];

	      total = 0;
	      for(int i=0; i< nval[S->ind[vari]]; i++)
		for(int j=0; j< nval[S->ind[varj]]; j++)
		  total += N[indx_start + i*step[vari] + j*step[varj]];

	      if(total > 0) {

		for(int i=0; i< nval[S->ind[vari]]; i++)
		  for(int j=0; j< nval[S->ind[varj]]; j++)
		    pairwise[i][j] = N[indx_start + i*step[vari] + j*step[varj]] / total;

		// Compute mutual information (gamma, the test statistic) for this distribution, using natural log
		for(int i=0; i< nval[S->ind[vari]]; i++)
		  marginal_vari[i] = 0;
		for(int j=0; j< nval[S->ind[varj]]; j++)
		  marginal_varj[j] = 0;

		for(int i=0; i< nval[S->ind[vari]]; i++) {
		  for(int j=0; j< nval[S->ind[varj]]; j++) {
		    marginal_vari[i] += pairwise[i][j];
		    marginal_varj[j] += pairwise[i][j];
		  }
		}

		mi = 0;
		for(int i=0; i< nval[S->ind[vari]]; i++) {
		  for(int j=0; j< nval[S->ind[varj]]; j++) {
		    if(pairwise[i][j] > 0) {
		      mi += pairwise[i][j] * log( pairwise[i][j] / (marginal_vari[i]*marginal_varj[j]) );
		    }
		  }
		}

		// Add (mi, total) to vector.
		MI_N_pairs.push_back(std::pair<double, int>(mi, total));
	      }

	      // increment conditioning set
	      done = true;
	      for (int i=0; i<S->n; i++) {
		if(i != vari && i != varj) {
		  if(assignment[i] < nval[S->ind[i]]-1) {
		    assignment[i]++;
		    done = false;
		    break;
		  }
		  else {
		    assignment[i] = 0;
		  }
		}
	      }
	    }

	    // Sort by mutual information in descending order
	    std::sort(MI_N_pairs.begin(), MI_N_pairs.end(), std::greater<pair<double,int> >());

	    double max_score_for_subset = 1; // initial value
	    double current_best = min_edge_scores[S->ind[vari]][S->ind[varj]];

	    // Iterate through the pairs and compute Beta
	    for(vector <std::pair <double,int> >::iterator it = MI_N_pairs.begin(); it != MI_N_pairs.end(); ++it) {

	      mi = (*it).first;
	      total = (*it).second;

	      // TEMP: if MI is greater than eta, just set beta to 1/2.
	      double log_beta = log_of_half;
	      if(mi < eta) {
	      
		// Compute beta -- call python  //TODO: integrate new C++ interpolation method right here!
		//pValue = PyObject_CallMethod(pInstance, "interpolateAndReturnLogBeta", "(idd)", total, mi, current_best);
		//Python method; interpolateAndReturnLogBeta(self, N, gamma, currSmallest=0)
		//if(PyErr_Occurred()) {
		//  cout << "Python error occured..." << total << "," << mi << " " << endl;
		//  PyErr_Print();
		//}
		//log_beta = PyFloat_AsDouble(pValue);
		
		//NEW: integration of C++ interpolation method
		log_beta = etasLgBetaTable.interpolate( total, mi, current_best); //double interpolate( int N, double gamma, double currSmallest = 0 )  
	      }

	      //		total_beta_calls++;
	      //		if(current_best == log_beta) stopped_early_beta_calls++;
	      //		if(total_beta_calls % 2 == 0) {
	      //		  double fract = stopped_early_beta_calls / total_beta_calls;
	      //		  cout << "Fraction of beta computations stopped early: " << fract << endl;
	      //		}

	      // Print!
	      //		cout << total << "," << mi << " " << beta << endl;

	      if(log_beta > max_score_for_subset || max_score_for_subset == 1) {
		max_score_for_subset = log_beta;

		// Quit early if possible
		if(max_score_for_subset >= current_best) {
		  break;
		}
		//		  else if(S->n >= 5) {
		//		    cout << "not quitting early: " << current_best << " " << max_score_for_subset << "  " << mi << "  " << total << endl;
		//		  }
	      }
	    }

	    //	    out << endl;
	    if(max_score_for_subset < current_best)
	      min_edge_scores[S->ind[vari]][S->ind[varj]] = max_score_for_subset;
	  }
	}
}

// ---------------------------------------------------
void model::learn_model(options &opt) {

	clock_t start_t = clock(); 
	
	double alpha = 1.0; opt.read("-alpha",alpha);
	int bic = opt.exists("-bic"); 
	int maxpa = 4; opt.read("-maxpa",maxpa); 
	maxpa = min(maxpa,nodes-1); 	
	int maxsize = maxpa+1; // largest subset we need to consider
	
	subset *S[maxsize+1]; 
	for (int n=0;n<=maxsize;n++) S[n] = new subset(nodes,n); 

	subset_array *SA[maxpa+1]; 
	SA[0] = new subset_array(nodes,0); 
	if (bic) SA[0]->set(S[0],bic_score(S[0])); 
	else SA[0]->set(S[0],score(S[0],alpha)); 
	
	VAR = new variable*[nodes]; 
	for (int i=0;i<nodes;i++) { VAR[i] = new variable; VAR[i]->id = i; }

	// Initialize edge score array. Will only populate in upper triangle. Initialize to 0.
	// When finished, will store min_S max_s ln[beta(n,gamma)]
	double** min_edge_scores;
	min_edge_scores = new double*[nodes];
	for(int i=0;i<nodes;i++) {
	  min_edge_scores[i] = new double[nodes];
	  for(int j=0;j<nodes;j++) {
	    min_edge_scores[i][j] = 0;
	  }
	}

	if(opt.exists("-edge_scores")) {

	  opt.read("-edge_scores",eta); assert(eta>0);
	  //	  char fname[LEN]; 
	  //	  opt.read("-edge_scores",fname); 
	  //	  ofstream out(fname);

	  // Load Python class for Beta computation
	  //PyObject *pName, *pModule, *pDict, *pClass, *pInstance, *pArgs;
	  //Py_Initialize(); // Initialize the Python Interpreter
	  //pName = PyString_FromString("betaTable"); // Build the name object
	  //pModule = PyImport_Import(pName); // Load the module object
	  //pDict = PyModule_GetDict(pModule); 
	  //pClass = PyDict_GetItemString(pDict, "betaTable"); // Build the name of a callable class
	  //pArgs = PyTuple_New(2);

	  // NOTE: if eta does not agree with the filename, will get bogus results. Have sanity check!
	  //PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble(eta));
	  string infilePath = "../interpolationData/";
	  string infileBase;
	  if(eta == .04) {
		  infileBase = infilePath + "etaFourHundredths"; 
		  //PyTuple_SetItem(pArgs, 1, PyString_FromString("/Users/eliotpbrenner/Documents/sontag/monteCarloBeta/test/etaFourHundredthsAdaptiveScalingLargeIncludingLargeGamma"));
	  }
	  else if(eta == .02) {
		infileBase = infilePath + "etaTwoHundredths"; 
	    //PyTuple_SetItem(pArgs, 1, PyString_FromString("/Users/eliotpbrenner/Documents/sontag/monteCarloBeta/test/etaTwoHundredthsAdaptiveScalingLargeIncludingLargeGamma"));
	  }
	  else if(eta == .01) {
		infileBase = infilePath + "etaOneHundredth"; 
	    //PyTuple_SetItem(pArgs, 1, PyString_FromString("/Users/eliotpbrenner/Documents/sontag/monteCarloBeta/test/etaOneHundredthAdaptiveScalingLargeIncludingLargeGamma"));
	  }
	  else if(eta == .005) {
		infileBase = infilePath + "etaFiveThousandths"; 
	    //PyTuple_SetItem(pArgs, 1, PyString_FromString("/Users/eliotpbrenner/Documents/sontag/monteCarloBeta/test/etaFiveThousandthsAdaptiveScalingLargeIncludingLargeGamma"));
	  }
	  else if(eta == .001) {
		infileBase = infilePath + "etaOneThousandth"; 
	//    PyTuple_SetItem(pArgs, 1, PyString_FromString("/Users/eliotpbrenner/Documents/sontag/monteCarloBeta/test/etaOneThousandthAdaptiveScalingLargeIncludingLargeGamma"));
	  }
	  else {
	    cout << "Invalid setting of eta." << endl;
	    return;
	  }

      //NEW: C++ interpolation method integration
	  string KL_div_gammaLTeta_InFile = infileBase + "N_KLDivPtsGammaLTEta.csv";
	  string lgBeta_gammaLTeta_InFile = infileBase + "LgBetaValsGammaLTEta.csv";
	  string KL_div_gammaGTeta_InFile = infileBase + "N_KLDivPtsGammaGTEta.csv";
	  string lgBeta_gammaGTeta_InFile = infileBase + "LgBetaValsGammaGTEta.csv";
	  //etasLgBetaTable is the persistent object for interpolation of lgBetas
	  LgBetaTable etasLgBetaTable(KL_div_gammaLTeta_InFile, lgBeta_gammaLTeta_InFile, 
						KL_div_gammaGTeta_InFile, lgBeta_gammaGTeta_InFile);

	  //pInstance = PyObject_CallObject(pClass, pArgs);
	  //PyErr_Print();
	  //cout << "loaded " << endl;
	  
	  for (int n=2;n<=min(MAX_CONDITIONING_SET_SIZE+2, min(maxsize+1, nodes));n++) { 
	    cout << "subsets size = " << n << endl;
	    subset *Se = new subset(nodes,n);
	    int cont = 1; Se->set(); 
	    while (cont) { 
	      edge_scores(Se, min_edge_scores, /*pInstance,*/ etasLgBetaTable);
	      cont = Se->next();
	    }
	  }

	  // Clean up Python
	  //Py_DECREF(pModule);
	  //Py_DECREF(pName);
	  //Py_Finalize();
	}

	clock_t end_t = clock(); 
	cout << "elapsed time = " << (end_t-start_t)/(double)CLOCKS_PER_SEC << " seconds so far. Now computing parent sets..." << endl;

	// Compute scores for the parent sets
	for (int n=1;n<=maxsize;n++) { 
		cout << "subsets size = " << n << endl; 
		if (n<maxsize) SA[n] = new subset_array(nodes,n); 
		int cont = 1; S[n]->set(); 
		while (cont) { 
			float val; if (bic) val = bic_score(S[n]); else val = score(S[n],alpha); 
			if (n<maxsize) SA[n]->set(S[n],val); 
			for (int i=0;i<n;i++) { 
				int id = S[n]->ind[i]; 
				S[n-1]->set(S[n],i); // set parents (all but i'th index) 
				// This next line subtracts off the parent's entropy from the set's entropy to obtain the
				// conditional entropy, which is precisely what is needed for the likelihood term. It also
				// correctly deals with the complexity penalty.
				float bscore = val - SA[n-1]->value(S[n-1]); 
				// parent sets for each variable are ordered in the increasing order of size
				parentset *tmp = new parentset(n-1,S[n-1]->ind,bscore); 
				tmp->prev = VAR[id]->paset; VAR[id]->paset = tmp;
			}			
			cont = S[n]->next(); 
		}
	}

	// Iterate over parent sets and add on edge scores
	cout << "adding on edge scores..." << endl;
	for (int i=0;i<nodes;i++) { 
	  variable *v = VAR[i];
	  for (parentset *s = v->paset; s != NULL; ) {
	    parentset *ps = s; s=s->prev; 

	    for (int j=0; j < ps->npa; j++) {
	      // Make ordered small to large
	      int first=i; int second=ps->pa[j];
	      if(first > second) {
		first=ps->pa[j];
		second=i;
	      }

	      ps->W += PSI2*min_edge_scores[first][second];
	    }
	  }
	}

	cout << "pruning parent sets" << endl; 
	for (int i=0;i<nodes;i++) { 
		// reverse the order in the list so that the smallest sets are on top
		variable *v = VAR[i]; parentset *PS = NULL; 
		for (parentset *s = v->paset; s != NULL; ) { parentset *tmp = s; s = s->prev; tmp->prev = PS; PS = tmp; }
		v->paset = PS; 
		for (int n=0;n<maxpa;n++) SA[n]->set(LOG_OF_ZERO); // pruning scores
		PS = NULL; unsigned long nps = 0; 
		for (parentset *s = v->paset; s != NULL; ) {
			parentset *ps = s; s=s->prev; 
			float W = LOG_OF_ZERO; // competition = max over smaller subsets
			for (int j=0;j<ps->npa;j++) W = max(W,SA[ps->npa-1]->value(ps->pa,ps->npa,j));
			SA[ps->npa]->set(ps->pa,ps->npa,max(W,ps->W)); // pruning criterion for larger sets
			if (ps->W > W) { // better than competition, keep
				ps->prev = PS; PS = ps; nps++; 
			} else { delete ps; } // prune
		}
		v->paset = PS; 
		cout << i << " pasets = " << nps << endl; 
	}
	for (int n=0;n<maxpa;n++) delete SA[n]; // don't need these any more
	
	end_t = clock(); 
	cout << "elapsed time = " << (end_t-start_t)/(double)CLOCKS_PER_SEC << " seconds total" << endl;

	// cleanup
	for(int i=0;i<nodes;i++)
	  delete min_edge_scores[i];
	delete min_edge_scores;
}

// ---------------------------------------------------
// TODO: Add sanity check that number of nodes equals number of items on a line (!)
void model::read_data(options &opt) {
	nodes = 0; opt.read("-nodes",nodes); assert(nodes>0);
	nval = new int[nodes]; for (int i=0;i<nodes;i++) nval[i] = 0; 
	
	// read data vectors
	char datafile[LEN]; opt.read("-data",datafile); 
	ifstream in(datafile); 
	D = NULL; 
	while (!in.eof()) {  
		data *d = new data(nodes); 
		d->read(in); in >> ws; // ws extracts whitespace characters.
		for (int i=0;i<nodes;i++) if (nval[i]<d->x[i]) nval[i] = d->x[i]; 
		d->prev = D; D = d; 
	}
	for (int i=0;i<nodes;i++) nval[i]++; 
}

// ---------------------------------------------------
void model::print_model(options &opt) {
	if (opt.exists("-mod_out")) {
		char fname[LEN]; 
		opt.read("-mod_out",fname); 
		ofstream out(fname);
		out.precision(dbl::digits10); // Set as high precision as necessary...
		out << nodes << endl; 
		for (int i=0;i<nodes;i++) VAR[i]->print(out); 
	}
}

// ---------------------------------------------------
model::model() { 
	nodes = 0; 
	nval = NULL;
	VAR = NULL; 
	D = NULL; 
	eta = 0;
}

// ---------------------------------------------------
model::~model() { 
	if (nodes>0) { 
		delete nval; 
	}
}


// ---------------------------------------------------
void model::calculate_epsilon(options &opt) {

  clock_t start_t = clock(); 
	
  int maxpa = 4; opt.read("-maxpa",maxpa); 
  maxpa = min(maxpa,nodes-1); 	
  int maxsize = maxpa+1; // largest subset we need to consider

  // Initialize edge score array. Will only populate in upper triangle. Initialize to 0.
  // When finished, will store min_S min_s tau
  double** min_edge_scores;
  min_edge_scores = new double*[nodes];
  for(int i=0;i<nodes;i++) {
    min_edge_scores[i] = new double[nodes];
    for(int j=0;j<nodes;j++) {
      min_edge_scores[i][j] = INFTY;
    }
  }


  for (int n=2;n<=min(MAX_CONDITIONING_SET_SIZE+2, min(maxsize+1, nodes));n++) { 
    cout << "subsets size = " << n << endl;
    subset *Se = new subset(nodes,n);
    int cont = 1; Se->set(); 
    while (cont) { 
      edge_taus(Se, min_edge_scores);
      cont = Se->next();
    }
  }

  // Output!
  char fname[LEN]; 
  opt.read("-epsilon",fname);
  ofstream out(fname);
  //  out << nodes << endl;

  for(int i=0;i<nodes;i++)
    for(int j=i+1;j<nodes;j++)
      out << i << " " << j << " " << min_edge_scores[i][j] << endl;

  clock_t end_t = clock();
  cout << "elapsed time = " << (end_t-start_t)/(double)CLOCKS_PER_SEC << " seconds." << endl;

  // cleanup
  for(int i=0;i<nodes;i++)
    delete min_edge_scores[i];
  delete min_edge_scores;
}


// ---------------------------------------------------
void model::edge_taus(subset* S, double** min_edge_scores) { 
	
	// prepare to index configurations
	int step[S->n+1]; int maxindx = 1; 
	for (int i=0;i<S->n;i++) { step[i]=maxindx; maxindx *= nval[S->ind[i]]; }
	
	double N[maxindx]; for (int i=0;i<maxindx;i++) N[i] = 0; 
	
	// go through all of the data, computing counts for this subset
	int ndata = 0;
	for (data *d=D;d != NULL;d=d->prev) { 
		int xindx = 0; for (int i=0;i<S->n;i++) xindx += d->x[S->ind[i]]*step[i]; 
		N[xindx] += 1.0; ndata++; 
	}

	// for all pairs i & j and for each setting of the other variables, compute tau
	int assignment[S->n];
	for (int vari=0; vari < S->n; vari++) {
	  for (int varj= vari+1; varj < S->n; varj++) {

	    double pairwise[nval[S->ind[vari]]][nval[S->ind[varj]]];
	    double marginal_vari[nval[S->ind[vari]]];
	    double marginal_varj[nval[S->ind[varj]]];

	    // initialize assignment to all zeros
	    for (int i=0; i<S->n; i++) assignment[i] = 0;

	    // Keep track of max tau for the subset
	    double max_tau_for_subset = -INFTY;

	    // for each value of the conditioning set
	    bool done = false;
	    while(!done) {

	      // compute the pairwise distribution
	      int indx_start = 0;
	      for (int i=0;i<S->n;i++)
		if(i!= vari && i != varj)
		  indx_start += assignment[i]*step[i];

	      double total = 0;
	      for(int i=0; i< nval[S->ind[vari]]; i++)
		for(int j=0; j< nval[S->ind[varj]]; j++)
		  total += N[indx_start + i*step[vari] + j*step[varj]];

	      if(total > 0) {

		for(int i=0; i< nval[S->ind[vari]]; i++)
		  for(int j=0; j< nval[S->ind[varj]]; j++)
		    pairwise[i][j] = N[indx_start + i*step[vari] + j*step[varj]] / total;

		// Compute mutual information (gamma, the test statistic) for this distribution, using natural log
		for(int i=0; i< nval[S->ind[vari]]; i++)
		  marginal_vari[i] = 0;
		for(int j=0; j< nval[S->ind[varj]]; j++)
		  marginal_varj[j] = 0;

		for(int i=0; i< nval[S->ind[vari]]; i++) {
		  for(int j=0; j< nval[S->ind[varj]]; j++) {
		    marginal_vari[i] += pairwise[i][j];
		    marginal_varj[j] += pairwise[i][j];
		  }
		}

		double mi = 0;
		for(int i=0; i< nval[S->ind[vari]]; i++) {
		  for(int j=0; j< nval[S->ind[varj]]; j++) {
		    if(pairwise[i][j] > 0) {
		      mi += pairwise[i][j] * log( pairwise[i][j] / (marginal_vari[i]*marginal_varj[j]) );
		    }
		  }
		}

		// Add (mi, total) to vector.
		if(mi > max_tau_for_subset)
		  max_tau_for_subset = mi;
	      }

	      // increment conditioning set
	      done = true;
	      for (int i=0; i<S->n; i++) {
		if(i != vari && i != varj) {
		  if(assignment[i] < nval[S->ind[i]]-1) {
		    assignment[i]++;
		    done = false;
		    break;
		  }
		  else {
		    assignment[i] = 0;
		  }
		}
	      }
	    }

	    if(max_tau_for_subset < min_edge_scores[S->ind[vari]][S->ind[varj]])
	      min_edge_scores[S->ind[vari]][S->ind[varj]] = max_tau_for_subset;

	  }
	}

}
