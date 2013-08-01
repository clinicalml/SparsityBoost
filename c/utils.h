#ifndef _UTILS_H
#define _UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cassert>

using namespace std;

// -----------------------------------

#include "options.h"

#ifndef NULL
#define NULL 0
#endif

// -----------------------------------
// elementary math utils

const double LOG_OF_ZERO = -1e+100;
const double MINUS_INFTY = -1e+100;
const double INFTY = 1e+100;

#define MAX(x,y) ( (x)>(y) ? (x) : (y) )
#define MIN(x,y) ( (x)<(y) ? (x) : (y) )

double maximum(double *a, int n); 
double sigmoid(double w0,double w1); 

// template <class T>
//inline T max(T a, T b) { 
//	if (a>b) return a; else return b; 
//}

//template <class T>
//inline T min(T a, T b) { 
//	if (a<b) return a; else return b; 
//}

// -----------------------------------
// matrix templates

template <class T>
T** new_matrix(int n, int m) { 
	assert(n>0 && m>0); 
	T** a = new T*[n]; 
	for (int i=0;i<n;i++) a[i] = new T[m]; 
	return a; 
}

template <class T>
void delete_matrix(T** a, int n, int m) { 
	for (int i=0;i<n;i++) delete a[i];
	delete a; 
}

template <class T>
void set_matrix(T** a, int n, int m, T val) { 
	for (int i=0;i<n;i++) 
		for (int j=0;j<m;j++) 
			a[i][j] = val;
}

template <class T>
void set_matrix(T** a, int n, int *nval, T val) { 
	for (int i=0;i<n;i++) 
		for (int j=0;j<nval[i];j++) 
			a[i][j] = val;
}

template <class T>
T** new_matrix(int n, int *nval) { 
	assert(n>0); 
	T** a = new T*[n]; 
	for (int i=0;i<n;i++) { a[i] = NULL; if (nval[i]>0) a[i] = new T[nval[i]]; }
	return a; 
}

template <class T>
void delete_matrix(T** a, int n, int *nval) { 
	for (int i=0;i<n;i++) if (nval[i]>0) delete a[i];
	delete a; 
}

template <class T>
void print(ostream& out, T** a, int n, int m) { 
	for (int i=0;i<n;i++) {
		for (int j=0;j<m-1;j++) out << a[i][j] << " "; 
		out << a[i][m-1] << endl;
	}
}

template <class T>
void print(ostream& out, T** a, int n, int *nval) { 
	for (int i=0;i<n;i++) {
		for (int j=0;j<nval[i]-1;j++) out << a[i][j] << " "; 
		if (nval[i]>0) out << a[i][nval[i]-1]; else out << "NULL";
		out << endl;
	}
}

template <class T>
void print(ostream& out, T* v, int n) { 
	for (int j=0;j<n-1;j++) out << v[j] << " "; 
	out << v[n-1] << endl;
}

// -----------------------------------
// random numbers

inline void randomize() { srand48( time(NULL) ); }
inline double randU() { return drand48(); }
inline int randI(int m, int n) 
{ return ( m+(int)((n-m+1)*randU()) % (n-m+1) ); }

// -----------------------------------
// sorting

// ----------------------------------------------
template <class T>
int sort_recur(T *A, int *ind, int *ind2, int ni, int nf) {
	
	if (ni>=nf) return 1; // nothing to do
	
	int k = ni+(nf-ni)/2; 
	sort_recur(A,ind,ind2,ni,k); 
	sort_recur(A,ind,ind2,k+1,nf); 
	
	int ki = ni; // counter for first half
	int kf = k+1; // counter for the second half
	
	// sort in ascending order to ind2
	for (int i=ni;i<=nf;i++)  
		if (ki>k) { ind2[i-ni] = ind[kf]; kf++; } // have to read from the second half
		else if (kf>nf) { ind2[i-ni] = ind[ki]; ki++; } // have to read from the first half
		else if (A[ind[ki]]>A[ind[kf]]) { ind2[i-ni] = ind[kf]; kf++; } // append smaller 
		else { ind2[i-ni] = ind[ki]; ki++; }
	
	// write back to ind
	for (int i=ni;i<=nf;i++) ind[i] = ind2[i-ni]; 
	
	return 1;
}

template <class T>
void sort(T *A, int *ind, int ni, int nf, int direc) {
	int *ind2 = new int[nf-ni+1]; 
	sort_recur(A,ind,ind2,ni,nf); // ascending 
	if (direc<0) { // descending
		for (int i=ni;i<=nf;i++) ind2[i-ni] = ind[i]; 
		for (int i=ni;i<=nf;i++) ind[i] = ind2[nf-i];
	}
	delete ind2;
} 

template <class T>
void Swap(T& a, T& b) { 
	T tmp = a; a = b; b = tmp;
}

#endif
