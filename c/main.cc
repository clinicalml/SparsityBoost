#include "utils.h"
#include "model.h"

#define LEN 4096

// ---------------------------------------------

int main(int argc, char* argv[]) {

	options opt(argc,argv);

	model *mod = new model();
	if (opt.exists("-data")) mod->read_data(opt);
	if (mod->nodes>0 && opt.exists("-epsilon")) {
	  mod->calculate_epsilon(opt); 
	  return 1;
	}
	if (mod->nodes>0) mod->learn_model(opt); 
	if (opt.exists("-mod_out") && mod->nodes>0) mod->print_model(opt); 
}
