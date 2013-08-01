#ifndef _LgBetaTable_H
#define _LgBetaTable_H

#include "Tau.h"
#include "binSearchObj.h"
#include "pwl_interp_2d.h"
#include <vector>

class LgBetaTable;
class LgBetaTable {
	protected:
		double eta;
		double t_eta;
		double p_eta[4];
		std::vector<int> N_int_LTEta;
		std::vector<int> N_int_GTEta;
		std::vector<int>::size_type numN_LTEta; 
		std::vector<int>::size_type numN_GTEta; 
		std::vector<std::vector<double> > KL_values_gammaLTeta;
		std::vector<std::vector<double> > lgBeta_values_gammaLTeta;
		std::vector<std::vector<double> > KL_values_gammaGTeta;
		std::vector<std::vector<double> > lgBeta_values_gammaGTeta;
		Tau TauCalc;  //TODO: make this constant and give binSearchObj a constructor which takes a const parameter
		const binSearchObj binSearchTau;
		double KL_from_p_eta(double t); 	
	public:
		LgBetaTable(const std::string &KL_div_gammaLTeta_InFile, 
				    const std::string &lgBeta_gammaLTeta_InFile, 
				    const std::string &KL_div_gammaGTeta_InFile, 
				    const std::string &lgBeta_gammaGTeta_InFile);
       double interpolate( int N, double gamma, double currSmallest = 0 );
	   double interpolateWithKL( int N, double KL, double gamma, double currSmallest = 0); //TODO: make protected
};
#endif
