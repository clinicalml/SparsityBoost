#ifndef _OPTIONS_H
#define _OPTIONS_H

#ifndef NULL 
#define NULL 0
#endif

class options_list;
class options_list {
public:
	char *opt; 
	options_list *next; 
	options_list() { opt = NULL; next = NULL; }
	~options_list();
};

class options { 
public: 
	options_list *OL; 
	int exists(char*);
	int read(char*,int&);
	int read(char*,double&);
	int read(char*,char*);
	void print(ostream &out);
	
	options() { OL = NULL; }
	options(int argc, char* argv[]);  
	options(istream& in); 
	~options();
};

#endif
