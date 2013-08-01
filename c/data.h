#ifndef _data_H
#define _data_H

#include "utils.h"

class data;
class data {
public:
	data *prev;
	int nvar;
	int *x;

	int read(istream& in);
	
	data(); 
	data(int nvar);
	~data();
};


#endif