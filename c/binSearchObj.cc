#include "binSearchObj.h"
#include "IncreasingFunctOnInterval.h"

#include <cmath>        // std::abs
binSearchObj::binSearchObj(double l, double u, IncreasingFunctOnInterval& f, double tol, int maxLevel) : lower(l), upper(u), f(f), tol(tol), maxLevel(maxLevel) {}

//main constructor to use: looks at the IncreasingFuctOnInterval to obtain its l and u bounds
binSearchObj::binSearchObj(IncreasingFunctOnInterval& f, double tol, int maxLevel) : lower(f.getLowerBound()), upper(f.getUpperBound()), f(f), tol(tol), maxLevel(maxLevel) {} 

double binSearchObj::search(double l, double u, double v, int currLevel) const {
	//printf("Calling Binary Search with l=%f, u=%f, v=%f, tol=%f\n",l,u,v,tol);
	if (currLevel > maxLevel) {
		//printf ("Returning because currLevel %d > maxlevel %d\n", currLevel, maxLevel );
		return l;
	}
	double m = (l+u)/2.0;
	double mVal = f.value(m);
	if (std::abs(mVal - v) < tol)  {
		//printf ("Returning because %f -%f smaller than %f\n", mVal, v,tol );
		return m;
	}
	else if (mVal > v) u = m;
	else if (mVal < v) l = m;
	return search(l,u,v,currLevel+1);
}
double binSearchObj::search(double v) const {
  double l = lower;
  double u = upper;
  return search(l, u, v, 0);
}
