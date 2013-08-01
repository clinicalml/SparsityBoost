#ifndef _Tau_H
#define _Tau_H
#include "IncreasingFunctOnInterval.h"

class Tau: public IncreasingFunctOnInterval{
	public:
		Tau();
		double value(double t) const;

};
#endif
