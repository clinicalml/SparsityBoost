#include "IncreasingFunctOnInterval.h"

void IncreasingFunctOnInterval::setBounds (float lb, float ub)
      { lowerBound=lb; upperBound=ub; }
	  IncreasingFunctOnInterval::IncreasingFunctOnInterval() {} //default constructor
	  IncreasingFunctOnInterval::~IncreasingFunctOnInterval() {}
double IncreasingFunctOnInterval::getLowerBound() { return lowerBound;}
double IncreasingFunctOnInterval::getUpperBound() { return upperBound;}
