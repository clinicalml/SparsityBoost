#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "LgBetaTable.h"
using namespace std;

typedef vector <double> record_t; //TODO: use these typedefs to make code below more readable
typedef vector <record_t> data_t;

//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.
istream& operator >> ( istream& ins, record_t& record )
  {
  // make sure that the returned record contains only the stuff we read now
  record.clear();

  // read the entire line into a string (a CSV record is terminated by a newline)
  string line;
  getline( ins, line );

  // now we'll use a stringstream to separate the fields out of the line
  stringstream ss( line );
  string field;
  while (getline( ss, field, ',' ))
    {
    // for each field we wish to convert it to a double
    // (since we require that the CSV contains nothing but floating-point values)
    stringstream fs( field );
    double f = 0.0;  // (default value is 0.0)
    fs >> f;

    // add the newly-converted field to the end of the record
    record.push_back( f );
    }

  // Now we have read a single line, converted into a list of fields, converted the fields
  // from strings to doubles, and stored the results in the argument record, so
  // we just return the argument stream as required for this kind of input overload function.
  return ins;
  }

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
istream& operator >> ( istream& ins, data_t& data )
  {
  // make sure that the returned data only contains the CSV data we read here
  data.clear();

  // For every record we can read from the file, append it to our resulting data
  record_t record;
  while (ins >> record)
    {
    data.push_back( record );
    }

  // Again, return the argument stream as required for this kind of input stream overload.
  return ins;  
  }

double interpolate1d(double x0,double x1,double y0,double y1,double x) {
	double slope = (y1-y0)/(x1-x0);
	return y0 + slope*(x-x0);
}

int readCSV(const string &CSVInfile, vector<std::vector<double> > &dataVector, vector<int> &N_int, double &eta, vector<int>::size_type &numN ) 
{ 
	/*
	 * Put the vector of data from the CSVInfile's 3rd line and onward into dataVector 
	 */

  data_t data;

  // Here is the file containing the data. Read it into data.
     
  ifstream infile(CSVInfile.c_str()); 
  infile >> data;

  // Complain if something went wrong.
  if (!infile.eof())
    {
    cout << "Read CSV Error " << CSVInfile << "!" << endl;
    return 1;
    }

  eta = data[0][0];

  
  N_int = std::vector<int>(data[1].begin(), data[1].end());
  //std::vector<int> N_int1(data[1].begin(), data[1].end());
  //N_int = N_int1;

  numN = N_int.size();
  int numNsRemoved = 0;
  for (unsigned n=0; n < numN; n++) {
	  vector<double> KL_row(data[n+2].begin() + 1, data[n+2].end());
	  if (KL_row.size() > 0)
	  {
	  dataVector.push_back(KL_row);
	  }
	  else
	  {
		  N_int.erase(N_int.begin()+n-numNsRemoved);
		  numNsRemoved++;
	  }
  }
  numN -= numNsRemoved;
  //for (unsigned n=0; n < N_int.size(); n++){
  //	  cout << N_int[n] << "\n";
  //}

  //cout << "Number of N is " << numN << endl;
  
  return 0;
}

double LgBetaTable::KL_from_p_eta(double t) {  //KL(p(t) || p_eta)  //TODO: for safety check for exceptional conditions
	double p_t[4] = {0.25 + t, 0.25 - t, 0.25 - t, 0.25 + t};
	return p_t[0]*log(p_t[0]/p_eta[0]) +  
		p_t[1]*log(p_t[1]/p_eta[1]) + 
		p_t[2]*log(p_t[2]/p_eta[2]) +
		p_t[3]*log(p_t[3]/p_eta[3]); 
		//TODO: could be optimized by combining 1st and 4th terms, and 2nd and 3rd terms
}
LgBetaTable::LgBetaTable(const string &KL_div_gammaLTeta_InFile, 
			const string &lgBeta_gammaLTeta_InFile, 
			const string &KL_div_gammaGTeta_InFile, 
			const string &lgBeta_gammaGTeta_InFile):
	binSearchTau(0.0,0.25,TauCalc, 0.0001,10)   //const data member initialized
{
	readCSV( KL_div_gammaLTeta_InFile, KL_values_gammaLTeta, N_int_LTEta, eta, numN_LTEta  );
	readCSV( lgBeta_gammaLTeta_InFile, lgBeta_values_gammaLTeta,N_int_LTEta, eta, numN_LTEta  );
	readCSV( KL_div_gammaGTeta_InFile, KL_values_gammaGTeta, N_int_GTEta,eta, numN_GTEta );
	readCSV( lgBeta_gammaGTeta_InFile, lgBeta_values_gammaGTeta,  N_int_GTEta, eta, numN_GTEta);
	t_eta = binSearchTau.search(eta); 
	p_eta[0] = 0.25 + t_eta; 
	p_eta[1] = 0.25 - t_eta;
	p_eta[2] = 0.25 - t_eta;
	p_eta[3] = 0.25 + t_eta;
}
		
double LgBetaTable::interpolate( int N, double gamma, double currSmallest /* = 0 */ ) 
{
   double t_gamma = binSearchTau.search(gamma);
   double currSmallestForCall = currSmallest;
   return min(interpolateWithKL(N,KL_from_p_eta(t_gamma),gamma, currSmallestForCall),currSmallest);  //middle arg is KL(p(t_gamma) || p_eta)
}
double LgBetaTable::interpolateWithKL( int N, double KL, double gamma, double currSmallest /* = 0 */ ) //TODO: make protected
{
	//cout << "Interpolating with N= " << N << ", KL= " << KL << ", gamma= " << gamma << endl;
  vector<int> *  N_int;
  vector<int>::size_type * numN;
  if (gamma > eta)  //TODO: adopt concise case notation with ?
  {
	  N_int = &(N_int_GTEta);
	  numN = &(numN_GTEta);
  }
  else
  {
	  N_int = &(N_int_LTEta);
	  numN = &(numN_LTEta);
  }
  //cout << "N_int->at(0): " << N_int->at(0) << endl;
  if (N < N_int->at(0)) {
	  return 0;  //given N is smaller than smallest tabulated N
  }

  int LgPos;  //position of first tabulated N greater than or equal to given N
  //test for case where N is larger than largest tabulated N
  //cout << "*numN - 1: " << *numN - 1 << endl;
  int lgTabN = N_int->at(*numN - 1);  //largest tabulated N
  //cout << "lgTabN: " << lgTabN << endl;
  if (N > lgTabN)
  {
	  LgPos = *numN - 1;
  }
  else //N is not larger than the largest tab N: do binary search for LgPos
  {
	  std::vector<int>::iterator N_low;
	  N_low = std::lower_bound (N_int->begin(), N_int->end(), N);
	  LgPos = N_low - N_int->begin();
  }
  int N_Lg = N_int->at(LgPos);
  //std::cout << "lower_bound at LgPos " << LgPos << "\n";

  vector< vector<double> > * p_KL_values;
  vector< vector<double> > * p_lgBeta_values;
  if (gamma > eta)
  {
	  p_KL_values = &(KL_values_gammaGTeta);
	  p_lgBeta_values = &(lgBeta_values_gammaGTeta);
  }
  else
  {
	  p_KL_values = &KL_values_gammaLTeta;
	  p_lgBeta_values = &(lgBeta_values_gammaLTeta);
  }

  if (N == N_Lg) //If N is actually in the list of tabulated N, only 1-d interpolation needed
  {
	  std::vector<double> KL_at_Pos = p_KL_values->at(LgPos);
	  std::vector<double> lgBeta_at_Pos = p_lgBeta_values->at(LgPos);
	  double KL_largest = KL_at_Pos.at(KL_at_Pos.size()-1);  //TODO: make function
	  if (KL > KL_largest) {
		  KL = KL_largest;
	  }
	  double smallest_KL = KL_at_Pos.at(0);
	  if (KL < smallest_KL) {
		  //cout << "KL less than smallest KL";
		  if (KL_at_Pos.size() < 2) return lgBeta_at_Pos.at(0);
		  double nextSmallest_KL = KL_at_Pos.at(1);
		  double first_LgBeta = lgBeta_at_Pos.at(0);
		  double second_LgBeta = lgBeta_at_Pos.at(1);
		  double * candidate = (gamma<eta)  ? &(second_LgBeta) : &(first_LgBeta);
		  if (currSmallest < *candidate) //TODO: eliminate the debug output and put this and similar if's on one line
		  {
		  //cout << "currSmallest " << currSmallest << " < " << *candidate << " , best candidate" << endl; 
		  return currSmallest;
		  }
		  return interpolate1d(smallest_KL, nextSmallest_KL,first_LgBeta,second_LgBeta,KL);
	  }
	  std::vector<double>::iterator KL_up;
	  KL_up = std::lower_bound (KL_at_Pos.begin(), KL_at_Pos.end(), KL);
	  int KL_LgPos = KL_up - KL_at_Pos.begin();
	  if (*KL_up == KL)  //unlikely, but have to test for it to avoid division by zero //TODO: absorb this check into the interpolate1d function
	  {
		  return lgBeta_at_Pos.at(KL_LgPos);
	  }
	  //1d interpolate
	  double x0 = KL_at_Pos.at(KL_LgPos-1); //TODO: add exception if KL_LgPos = 0
	  double x1 = *KL_up;
	  double x = KL;
	  double y0 = lgBeta_at_Pos.at(KL_LgPos-1); //TODO: add exception if KL_LgPos = 0
	  double y1 = lgBeta_at_Pos.at(KL_LgPos);
	  double * candidate = (gamma<eta)  ? &(y1) : &(y0);
	  if (currSmallest < *candidate) 
	  {
		  //cout << "currSmallest " << currSmallest << " < " << *candidate << " , best candidate" << endl; 
		  return currSmallest;
	  }
	  //cout << "currSmallest " << currSmallest << " > " << *candidate << ", best candidate" << endl; 
	  //cout << "1d Interpolation with " << " x0: " << x0 <<" x1: " << x1 << " y0: " << y0 <<" y1: " << y1  <<" x: " << x << endl;
	  return interpolate1d(x0,x1,y0,y1,x);
  
  }  //end of 1-d interpolation block to execute if given N is in list of tabulated N
  else {  //2-d interpolation needed
	  int SmPos = LgPos - 1 ;  //the position smaller than LgPos
	  if (SmPos < 0) //test for error situation: if we get to this point LgPos should be at least 1
	  {
		  return 0;
	  }
	  //get KL and LgBeta values at LgPos and SmPos
	  std::vector<double> KL_at_LgPos = p_KL_values->at(LgPos);
	  std::vector<double> lgBeta_at_LgPos = p_lgBeta_values->at(LgPos);
	  std::vector<double> KL_at_SmPos = p_KL_values->at(SmPos);
	  std::vector<double> lgBeta_at_SmPos = p_lgBeta_values->at(SmPos);

	  //Thresholding: if test_KL larger than largest KL stored for the smaller N, replace with this largest KL
	  double KL_largest = KL_at_SmPos[KL_at_SmPos.size()-1];  //TODO: make function
	  if (KL > KL_largest) {
		  KL = KL_largest;
	  }
	  //Find the values to interpolate 
	  std::vector<double>::iterator KL_up_SmN, KL_up_LgN; 
	  KL_up_SmN = std::lower_bound (KL_at_SmPos.begin(), KL_at_SmPos.end(), KL);
	  KL_up_LgN = std::lower_bound (KL_at_LgPos.begin(), KL_at_LgPos.end(), KL);
	  
	  //another form of thresholding: if these iterators point to beginning of the vectors, increment them
	  if (KL_up_SmN == KL_at_SmPos.begin()) 
	  {  //one, therefore both, iterators point to the start of the vectors: iterate them
		  ++KL_up_SmN;
		  ++KL_up_LgN;
	  }
	  int KL_LgPos_SmN = KL_up_SmN - KL_at_SmPos.begin(); 
	  int KL_LgPos_LgN = KL_up_LgN - KL_at_LgPos.begin(); 
	  double KL0_SmN = KL_at_SmPos.at(KL_LgPos_SmN-1); //TODO: add exception if KL_LgPos_SmN == 0
	  double KL1_SmN = *KL_up_SmN;
	  double KL0_LgN = KL_at_SmPos.at(KL_LgPos_LgN-1); //TODO: add exception if KL_LgPos_SmN == 0
	  double KL1_LgN = *KL_up_LgN;
	  double LgBeta0_SmN = lgBeta_at_SmPos[KL_LgPos_SmN-1]; 
	  double LgBeta1_SmN = lgBeta_at_SmPos[KL_LgPos_SmN];
	  double LgBeta0_LgN = lgBeta_at_LgPos[KL_LgPos_LgN-1];
	  double LgBeta1_LgN = lgBeta_at_LgPos[KL_LgPos_LgN];
	  double * candidate = (gamma<eta)  ? &(LgBeta1_LgN) : &(LgBeta0_LgN);
	  if (currSmallest < *candidate) 
	  {
		  //cout << "currSmallest " << currSmallest << " < " << *candidate << " , best candidate" << endl; 
		  return currSmallest;
	  }
	  //cout << "currSmallest " << currSmallest << " > " << *candidate << ", best candidate" << endl; 
	  int N_Sm = N_int->at(LgPos-1); //TODO: add exception if LgPos = 0
	  //cout << "Interpolating the data N = " << N << " KL " << KL << endl;
	  //cout << "Over the points  1     N = " << N_Sm << " KL " << KL0_SmN << " LgBeta " << LgBeta0_SmN << endl; 
	  //cout << "Over the points  2     N = " << N_Sm << " KL " << KL1_SmN << " LgBeta " << LgBeta1_SmN << endl; 
	  //cout << "Over the points  3     N = " << N_Lg << " KL " << KL0_LgN << " LgBeta " << LgBeta0_LgN << endl; 
	  //cout << "Over the points  4     N = " << N_Lg << " KL " << KL1_LgN << " LgBeta " << LgBeta1_LgN << endl; 
	  int nxd = 2; //number of x grid points
	  int nyd = 2; //  ...     y ...
	  int ni = 1; //number of interpolant points
	  double xd[2] = {N_Sm, N_Lg};
	  double yd[2] = {KL0_SmN, KL1_SmN};
	  double zd[4] = { LgBeta0_SmN, LgBeta0_LgN,LgBeta1_SmN,  LgBeta1_LgN }; 
	  double xi[1] = {N*1.0};
	  double yi[1] = {KL};
	  double *interpolatedLgBeta;
	  interpolatedLgBeta = pwl_interp_2d (nxd, nyd, xd, yd, zd, ni, xi, yi);
	  return *interpolatedLgBeta;  //NB: rely on the public calling function to take the minimimum with the currSmallest
  }

}


//-----------------------------------------------------------------------------
// Now to put it all to use.

/*  TODO: put this testing code into a proper unit testing file
int main(int argc, char* argv[]) 
{
	if (argc < 3) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
        std::cout << "Usage is -in <infileBase>\n"; // Inform the user of how to use the program
        std::cin.get();
        exit(0);
    } else 
	{ // if we got enough parameters...
        string infileBase;
        std::cout << argv[0];
        for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                                          * Note that we're starting on 1 because we don't need to know the 
                                          * path of the program, which is stored in argv[0] */
 /*               if (string(argv[i]) == "-in") 
				{
                    // We know the next argument *should* be the filename:
                    infileBase = string(argv[++i]);
                } 
                else 
				{
                    std::cout << "   " << i << argv[i] << "  Not enough or invalid arguments, please try again.\n";
					//Sleep(200); 
                    exit(0);
				}
            std::cout << argv[i] << " ";
        }
	  // Here is the data we want.
	  // TODO: generalize to gamma greater than eta
	  
	  std::vector<int> testNVect;
	  testNVect.push_back(100);
	  testNVect.push_back(200);
	  testNVect.push_back(99);
	  testNVect.push_back(999);
	  testNVect.push_back(8999);
	  testNVect.push_back(10000);
	  testNVect.push_back(11000);
	  testNVect.push_back(12000);
	  testNVect.push_back(13000);
	  testNVect.push_back(14000);
	  testNVect.push_back(15000); 
	  
	  std::vector<double> testGammaVect;
	  testGammaVect.push_back(1e-4);
	  testGammaVect.push_back(6e-4);
	  testGammaVect.push_back(1e-3);
	  testGammaVect.push_back(10e-3);
	  testGammaVect.push_back(100e-3);
	  testGammaVect.push_back(200e-3);
	  testGammaVect.push_back(300e-3);
	  
	  
	  string KL_div_gammaLTeta_InFile = infileBase + "N_KLDivPtsGammaLTEta.csv";
	  string lgBeta_gammaLTeta_InFile = infileBase + "LgBetaValsGammaLTEta.csv";
	  string KL_div_gammaGTeta_InFile = infileBase + "N_KLDivPtsGammaGTEta.csv";
	  string lgBeta_gammaGTeta_InFile = infileBase + "LgBetaValsGammaGTEta.csv";
	  LgBetaTable etasLgBetaTable(KL_div_gammaLTeta_InFile, lgBeta_gammaLTeta_InFile, 
						KL_div_gammaGTeta_InFile, lgBeta_gammaGTeta_InFile);


	  for (std::vector<int>::iterator it = testNVect.begin(); it != testNVect.end(); ++it) 
	  {
		  for (std::vector<double>::iterator it2 = testGammaVect.begin(); it2 != testGammaVect.end(); ++it2) {
			  cout << endl;
			  cout << "******************************************************************" << endl;
			  std::cout << "N: " << *it << ",  Gamma: " << *it2 << endl;
			  double interpolant =  etasLgBetaTable.interpolate( *it, *it2, -5 );
			  cout << "Interpolant: " << interpolant  << endl;
		  }
	  } 
	  
	   cout << "Good bye!\n";
	  return 0;
	}
}
*/
