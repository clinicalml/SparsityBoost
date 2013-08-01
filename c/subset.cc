#include "subset.h"

// ---------------------------------------------------
subset::subset(int nodes, int n) { 
	this->nodes = nodes; 
	this->n = n; 
	ind = new int[n+1]; 
	assert(n<=nodes); 
	for (int i=0;i<n;i++) ind[i] = i; // set to first n
	ind[n] = nodes; 
}

// ---------------------------------------------------
subset::~subset() { 
	if (n>0) delete ind; n = 0; nodes = 0; 
}

// ---------------------------------------------------
void subset::set() { 
	for (int i=0;i<n;i++) ind[i] = i; ind[n] = nodes; 
}

// ---------------------------------------------------
void subset::set(subset *S) { 
	assert(S->n == n && S->nodes == nodes); 
	for (int i=0;i<n;i++) ind[i] = S->ind[i]; ind[n] = nodes; 
}

// ---------------------------------------------------
void subset::set(subset *S, int k) { 
	assert(S->n == n+1 && S->nodes == nodes); 
	int s = 0; for (int i=0;i<S->n;i++) if (i != k) { ind[s++] = S->ind[i]; } 
	ind[n] = nodes; 
}

// ---------------------------------------------------
int subset::next() { 
	int ok = 0; 
	for (int s = n-1; !ok && s>=0; s--) { 
		if (ind[s]+1<ind[s+1]) {
			ind[s]++; // advance this position
			for (int i=s+1;i<n;i++) ind[i] = ind[s]+i-s; // re-initialize all higher indexes
			ok = 1; 
		}
	}
	return ok; // ok = 0 if we cannot advance further
}

// ---------------------------------------------------
// nodes choose n
unsigned long F(int nodes, int n) { 
	if (n<=0) return 1; 
	unsigned long index = 1; 
	for (int i=0;i<n;i++) index *= (nodes-i); 
	for (int i=0;i<n;i++) index /= (i+1); 
	return index; 
}

// ---------------------------------------------------
subset_array::subset_array(int nodes, int n) { 
	assert(nodes>0 && n<=nodes); 
	this->nodes = nodes; this->n = n; 
	
	if (n==0) { 
		maxindex = 1; 
		A = new float[maxindex]; A[0] = 0; 
		step = NULL;
	} else { 
		maxindex = F(nodes,n); 
		A = new float[maxindex]; // one element per subset
		for (unsigned long s=0; s<maxindex; s++) A[s] = 0; 
		
		step = new unsigned long*[n]; 
		for (int i=0;i<n;i++) { 
			step[i] = new unsigned long[nodes-n+1]; 
			step[i][0] = 0; 
			for (int j=1;j<nodes-n+1;j++) step[i][j] =  F(nodes-i-1-j,n-i-1)+step[i][j-1]; 
		}
//		unsigned long index = 0; for (int i=0;i<n;i++) index += step[i][nodes-n]; 
//		cout << "nodes = " << nodes << " n = " << n << " " << maxindex << " " << index << endl;
//		assert(index == maxindex-1); 
	}
}
// ---------------------------------------------------
subset_array::~subset_array() { 
	delete A; 
	if (step != NULL) { 
		for (int i=0;i<n;i++) delete step[i]; delete step; 
	}
}

// ---------------------------------------------------
float subset_array::value(subset *S) { 
	assert(S->n == n); 
	unsigned long index = 0; 
	for (int i=0;i<n;i++) index += step[i][S->ind[i]-i]; 
	return A[index]; 
}

// ---------------------------------------------------
float subset_array::value(int *I, int ni, int k) { 
	assert(ni == n+1); 
	unsigned long index = 0; 
	for (int i=0;i<n;i++) { 
		int s = i; if (i>=k) s = i+1; // skip k'th entry in I
		index += step[i][I[s]-i]; 
	}
	return A[index]; 
}

// ---------------------------------------------------
void subset_array::set(subset *S, float val) { 
	assert(S->n == n); 
	unsigned long index = 0; 
	for (int i=0;i<n;i++) index += step[i][S->ind[i]-i]; 
	A[index] = val; 
}

// ---------------------------------------------------
void subset_array::set(int *I, int ni, float val) { 
	assert(ni == n); 
	unsigned long index = 0; 
	for (int i=0;i<n;i++) index += step[i][I[i]-i]; 
	A[index] = val; 
}

// ---------------------------------------------------
void subset_array::set(float val) { 
	for (unsigned long s=0;s<maxindex;s++) A[s] = val; 
}
