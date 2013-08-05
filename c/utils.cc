#include "utils.h"

// ---------------------------------------------------
double maximum(double *a, int n) { 
	int ind = 0; 
	for (int i=1;i<n;i++) if (a[i]>a[ind]) ind = i; 
	return a[ind]; 
}

// ---------------------------------------------------
double sigmoid(double w0, double w1) { 
	if (w0>w1) { 
		double p = exp(w1-w0); 
		return p/(p+1.0); 
	} else { 
		double p = exp(w0-w1); 
		return 1.0/(p+1.0); 
	}
}
