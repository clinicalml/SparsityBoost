#ifndef _IncreasingFunctOnInterval_H
#define _IncreasingFunctOnInterval_H

class IncreasingFunctOnInterval;
class IncreasingFunctOnInterval {
  protected:
    double lowerBound, upperBound;
  public:
    IncreasingFunctOnInterval(); //default constructor
	~IncreasingFunctOnInterval();
	void setBounds (float lb, float ub);
    double getLowerBound();
	double getUpperBound();
	virtual double value (double input) const =0;
  };

#endif
