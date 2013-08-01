#include "Tau.h"
#include "IncreasingFunctOnInterval.h"
#include <math.h>       /* log */

Tau::Tau()
			: IncreasingFunctOnInterval( ) //Call default constr. of base
	{
		setBounds(0.0,0.24999);  //t ranges from 0 to 0.25=t_{\max}
	} //constructor has empty body

double Tau::value(double t) const
{ if (t< lowerBound)
	t = lowerBound;
if (t> upperBound)
	t = upperBound;
return (0.5+2*t)*log(1+4*t)+(0.5-2*t)*log(1-4*t);
}



