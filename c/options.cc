#include "utils.h"
#include "options.h"
#include <cstring>

#define STRLEN 4096

// --------------------------------------
options_list::~options_list() { 
	if (next != NULL) delete next; 
	next = NULL; 
	if (opt != NULL) delete opt; 
	opt = NULL;
}

// --------------------------------------
void options::print(ostream &out) {
	for (options_list *ol=OL; ol != NULL; ol=ol->next) {
		out << ol->opt << " "; 
	}
	out << endl;
}

// --------------------------------------
options::options(int nargc, char* argv[]) {
	OL = new options_list; 
	options_list *ol = OL;
	for (int i=0;i<nargc;i++) {
		if (ol->opt != NULL) { ol->next = new options_list; ol = ol->next; }
		ol->opt = new char[strlen(argv[i])+1];
		strcpy(ol->opt,argv[i]);
	}
}

// --------------------------------------
options::options(istream &in) { 
	char linebuf[STRLEN];
	in >> ws; in.getline(linebuf,STRLEN); 
	istringstream iss(linebuf); 
	
	OL = new options_list; 
	options_list *ol = OL;
	while (!iss.eof()) { 
		if (ol->opt != NULL) { ol->next = new options_list; ol = ol->next; }
		char buf[STRLEN]; iss >> buf >> ws; 	
		ol->opt = new char[strlen(buf)+1];
		strcpy(ol->opt,buf);
	}	
}

// --------------------------------------
int options::exists(char *id) {
	int found = 0;
	for (options_list *ol=OL; ol != NULL && found==0; ol=ol->next)
		if (strcmp(id,ol->opt)==0) found = 1; 

	return found;
}

// --------------------------------------
int options::read(char* id, int& n) {
	int found = 0;
	for (options_list *ol=OL; ol != NULL && found==0; ol=ol->next)
		if (strcmp(id,ol->opt)==0 && ol->next != NULL) { 
			sscanf(ol->next->opt,"%d",&n); 
			found = 1;
		}
		
	return found;
}

// --------------------------------------
int options::read(char* id, double& x) {
	int found = 0;
	for (options_list *ol=OL; ol != NULL && found==0; ol=ol->next)
		if (strcmp(id,ol->opt)==0 && ol->next != NULL) { 
			sscanf(ol->next->opt,"%lf",&x); 
			found = 1;
		}
		
	return found;
}

// --------------------------------------
int options::read(char* id, char *str) {
	int found = 0;
	for (options_list *ol=OL; ol != NULL && found==0; ol=ol->next)
		if (strcmp(id,ol->opt)==0 && ol->next != NULL) { 
			strcpy(str,ol->next->opt);
			found = 1;
		}
		
	return found;
}

// --------------------------------------

options::~options() {
	if (OL != NULL) delete OL; 
	OL = NULL;
}

