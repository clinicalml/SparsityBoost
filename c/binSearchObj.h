#ifndef _binSearchObj_H
#define _binSearchObj_H
#include "IncreasingFunctOnInterval.h"

class binSearchObj;
class binSearchObj {
  private:
	  const double lower,upper; //fixed upper and lower bounds for the binSearch
	  const double tol; //for terminating search
	  const int maxLevel; //maxlevel of recursion
	  const IncreasingFunctOnInterval& f; //because abstract class, must be ref.

  public:
  binSearchObj(double l, double u, IncreasingFunctOnInterval& f, double tol, int maxLevel);
  
  //main constructor to use: looks at the IncreasingFuctOnInterval to obtain its l and u bounds
  binSearchObj(IncreasingFunctOnInterval& f, double tol, int maxLevel);

  double search(double l, double u, double v, int currLevel) const;
  double search(double v) const;
};
#endif
