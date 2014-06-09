/**
 *  
 * @Description This program computes the profile matrix of a protein fragment with respect 
* to the twenty natural amino acid types. For each residue position in the 
* fragment, amino acid types are ordered by tap score. The profile matrix is 
* printed to 'tap_profile.txt'. A summary of the matrix is printed both at the 
* end of the matrix as well as on the screen. The summary is printed in the 
* form of a table, whose column "<=N" indicates the number of residue 
* positions for which the native amino acid type has a tap score in the 
* highest N (N=1,...,5). Column "N>5"  indicates the number of remaining 
* residue positions.
*/
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <Potential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmega.h>
#include <PhiPsiOmegaChi1.h>
#include <PhiPsiOmegaChi1Chi2.h>
#include <Omega.h>
#include <Chi1Chi2.h>
#include "Sort.h"

using namespace Biopool;


/* Type definition for a record that binds a tap score to an 
   amino acid type. Relational operators, based on the tap 
   score only, are also defined to allow a sorting of these 
   records (i.e., these operators satisfy the conditions 
   specified in Sort.h for the generic sortable item). */
struct TapRec {
  double tap;
  AminoAcidCode amino; // amino acid type

  TapRec() {}
  TapRec(const double t, AminoAcidCode a) : tap(t), amino(a) {}
};


bool operator==(const TapRec& t1, const TapRec& t2) {
  return t1.tap == t2.tap;  
}


bool operator<(const TapRec& t1, const TapRec& t2) {
  return t1.tap < t2.tap;  
}


bool operator>(const TapRec& t1, const TapRec& t2) {
  return t1.tap > t2.tap;  
}


///////////////////////////////////////////////////////////////


void sShowHelp(){
    cout << "tap4design - Protein design using TAP scores. \n\n"
       << "Developed by Andrea Bazzoli <bazzoli@dti.unimi.it> \n"
       << "\nOPTIONS:\n"
       << "\t-i <filename> \t\t Input PDB file." <<"\n"
       << "\t[--chain <id>]  \t ID of chain to load from PDB file\n"
       << "\t[-s <int> ] \t\t Start residue (DEFAULT = second AA)\n"
       << "\t[-e <int> ] \t\t End residue (DEFAULT = last-1)\n"
       << "\n\t--help  \t\t Extended help for the method\n"
       << endl;
}


void sShowSpecial(){
  sShowHelp();

  cout << "\nADDITIONAL OPTIONS:\n"
       << "\t[-v]    \t\t Verbose mode\n"
       << "\t[--know <filename>] \t Knowledge filename (Default = tor.par).\n"
       << "\t[--pp]      \t\t consider Phi and Psi Angles\n"
       << "\t[--ppo]     \t\t consider Phi, Psi, Omega Angles\n"
       << "\t[--ppoc]    \t\t consider Phi, Psi, Omega, and Chi1 Angles\n"
       << "\t[--ppocc]   \t\t consider Phi, Psi, Omega, Chi1, and Chi2 Angles (DEFAULT) \n"
       << "\t[--bin <integer>] \t the bin for Phi and Psi angles (Default = 10 degrees).\n"
       << endl;
}


/**
 *@Description returns the tap score
*@param       
* diheds: a set of backbone plus side-chain, if any, dihedral angles.
* code: an amino acid type.
*   pot: a torsional-potential function
*@return     
*    tap score, calculated with 'pot', of the amino acid whose dihedral 
*    angles are specified by 'diheds' and whose type is specified by code.
*/              
double getTap(AminoAcid& diheds, AminoAcidCode code, TorsionPotential* pot) {

    double eMin = -log(pot->pReturnMinPropensities(code));
    double eMax = -log(pot->pReturnMaxPropensities(code));

    // energy of the assembled amino acid according to 'pot'
    double energy = pot->calculateEnergy(diheds, code);

    return (energy - eMin) / (eMax - eMin);
}


/**
 *@Description Builds a vector of tap records, which are specific to different amino 
*    acid types but refer to the same configuration of dihedral angles.
*    The tap record of the i-th natural amino acid type in 'AminoAcidCode' is 
*    stored in the i-th position of the vector (i=0,..,19).
*@return   pointer to the vector of tap records.
 *@param    
*    diheds: an amino acid containing the fixed configuration
*            of dihedral angles.
*    pot: pointer to the torsional-potential function used to compute tap 
*         scores.
*
*  Notes: 
*    this function relies on the fact that the codes of the 20 natural amino
*    acid types are 20 consecutive integers.
*    the vector of tap records is allocated from the heap memory and should be
*    deleted explicitly by some code outside this function.
*/
vector<TapRec>* getAllAminoTaps(AminoAcid& diheds, TorsionPotential* pot) {
    vector<TapRec>* aminoTaps = new vector<TapRec>;
    aminoTaps->reserve(AminoAcid_CODE_SIZE-1);
    
    for(int j=ALA; j<=TYR; ++j) 
    {
        double tap = getTap(diheds, AminoAcidCode(j), pot);
        TapRec tr(tap, AminoAcidCode(j));

        (*aminoTaps).push_back(tr);
    }

    return aminoTaps;
}


/**
 *@Description  builds the tap profile matrix of a protein-chain fragment. The i-th column 
*    of this matrix contains the tap scores of the 20 natural amino acid types  
*    for the i-th position of the fragment (i=0,...,fragLen-1). For each 
*    column, the j-th row contains the tap record of the j-th natural amino  
*    acid type in the AminoAcidCode enumeration.
*@Return  pointer to the profile matrix.
*@param  chain: a Spacer object containing the entire protein chain.
* start: index, in the chain, of the fragment's start residue.
* end: index, in the chain, of the first residue after the fragment.
 * pot: pointer to the torsional-potential function used to compute tap  scores.
* Notes: 
*    the profile matrix is allocated as a vector V of pointers to vectors, so
*    the i-th column of the matrix is actually "pointed to" by the i-th 
*    element of V (i=0,...,fragLen-1).
*    the profile matrix is allocated from the heap memory and should be 
*    deleted explicitly by code outside this function.
*/
vector<vector<TapRec>*>* buildFragTapProfileMatrix(Spacer& chain, 
  const unsigned int start, const unsigned int end, TorsionPotential* pot) 
{
  vector<vector<TapRec>*>* profMat = new vector<vector<TapRec>*>;
  profMat->reserve(end - start);

  for(unsigned int i=start; i<end; ++i)
  {
    vector<TapRec>* col = getAllAminoTaps(chain.getAmino(i), pot); 
    profMat->push_back(col);
  }

  return profMat;
}


/**
 * @description sorts a vector of tap records.
*@param pointer to the vector.
* Notes: 
*    this function calls the SortAlgos::quicksort() template function to do
*    the sorting. Since the latter operates on arrays, the vector is first 
*    copied to an array and then the sorted array is copied back to the vector.
*/
void sortTapRecVec(vector<TapRec>* v) 
{
    const int size = v->size();
    TapRec* a = new TapRec[size];

    // copy vector to array
    for(int i=0; i<size; ++i)
      a[i] = (*v)[i];

    SortAlgos::quicksort(a, 0, size-1);

    // copy sorted array to vector
    for(int i=0; i<size; ++i)
      (*v)[i] = a[i];

    delete[] a;
}


/**
 * @description  prints to file the tap profile matrix of a protein-chain fragment. 
*    Amino acid types are printed using the three-letter code. Each amino acid
*    type is followed by its tap score.
*@param mtx: pointer to the profile matrix.
*      chain: a Spacer object containing the entire protein chain.
*      start: index, in the chain, of the fragment's start residue.
*      end: index, in the chain, of the first residue after the fragment.
*      fileName: name of the output file.
*      Notes: 
*       since a tap score equal to zero can be printed as -0.000, I chose to print
*      the absolute values of tap scores. (Remember that tap scores are always
*      non-negative.)
*      The j-th row of the matrix is printed as the j-th row in the file 
*      (j=0,...,19). For each row, the i-th tap record in the matrix is printed 
*      as the i-th tap record in the file. Therefore, the i-th column of the 
*      printed matrix describes the tap profile of the i-th residue in the 
*      fragment. A header row is printed, and its i-th item, an extension of the
*      i-th column of the matrix, also identifies the i-th residue in the 
*      fragment. (i=0,...,fragLen-1.)
*/
void printFragTapProfileMatrix3L(
  vector<vector<TapRec>*>* mtx, Spacer& chain, const unsigned int start, 
  const unsigned int end, const char* fileName)
{

  // open output file stream
  std::ofstream outfs(fileName);
  if(!outfs)
    ERROR("Can't open " + string(fileName) + ".", exception);

  // print header line
  outfs.fill('#');
  for(unsigned int i=start; i<end; ++i)
  {
    string iCode = aminoAcidThreeLetterTranslator(
                     AminoAcidCode(chain.getAmino(i).getCode()));
    outfs << iCode;
    outfs << setw(6);
    outfs << i << "  ";
  }
  outfs << "\n\n";

  // change printing format for floating-point values 
  unsigned int oldPrec = outfs.precision();
  ios::fmtflags oldFlags = outfs.flags();
  outfs.setf(ios::fixed, ios::floatfield);

  // print one row of the matrix at a time
  const int ROWS = AminoAcid_CODE_SIZE - 1;
  for(int j=0; j<ROWS; ++j)
  {
    const unsigned int COLS = end - start;
    for(unsigned int i=0; i<COLS; ++i)
    {
      AminoAcidCode iCode = (*(*mtx)[i])[j].amino; 
      double iTap = (*(*mtx)[i])[j].tap;

      // separator character between amino acid type and tap score. Use a special
      // character for the native amino acid type
      AminoAcidCode nCode = AminoAcidCode(chain.getAmino(start+i).getCode());
      char sep = (iCode == nCode) ? '#' : '_';

      outfs << aminoAcidThreeLetterTranslator(iCode) 
            << sep 
            << setw(5) << setprecision(3) << fabs(iTap) << "  ";      
    }

    outfs << "\n";
  }

  outfs << "\n\n";

  // reset the floating-point printing format to its original status
  outfs.precision(oldPrec);
  outfs.flags(oldFlags);
} 


/**
* @param   mtx: pointer to the profile matrix of a protein-chain fragment.
*          chain: a Spacer object containing the entire protein chain.
*          start: index, in the chain, of the fragment's start residue.
*          end: index, in the chain, of the first residue after the fragment.
*          row: index of a row in the matrix.
* @Return value:
*    number of columns (residue positions) having a native amino acid type in 
*    that row.
*/
int getRowNativeCoverage(vector<vector<TapRec>*>* mtx, Spacer& chain, 
  const unsigned int start, const unsigned int end, const int row) 
{
  int tot = 0;
  const int COLS = end - start;
  for(int i=0; i<COLS; ++i) {

    AminoAcidCode nCode = AminoAcidCode(chain.getAmino(start+i).getCode());
    AminoAcidCode iCode = (*(*mtx)[i])[row].amino;

    if(nCode == iCode)
      tot++; 
  }
  
  return tot;
}


/**
 * @Author:         Andrea Bazzoli
* @Description:
*    prints the summary of the tap-score profile matrix for a protein-chain fragment.  
*@param    
*    mtx: pointer to the profile matrix of the protein-chain fragment.
*    chain: a Spacer object containing the entire protein chain.
*     start: index, in the chain, of the fragment's start residue.
*    end: index, in the chain, of the first residue after the fragment.
*    fileName: name of the output file. If this argument is omitted, then the 
*              summary is printed on standard output.
* Notes: this function uses the non-portable identifier ios::app to open 
*         an output file stream in append mode.
*/
void printSummaryFTPM(vector<vector<TapRec>*>* mtx, Spacer& chain, 
  const unsigned int start, const unsigned int end, const char* fileName=0) 
{
  
  // select whether to direct output to file or screen
  ostream* tapOut = 0;
  if(fileName) 
  {
    tapOut = new std::ofstream(fileName, ios::app); 
    if(!(*tapOut))
      ERROR("Can't open " + string(fileName) + ".", exception);
  }
  else
    tapOut = &cout;

  // print header line
  (*tapOut).fill('#');
  for(int i=1; i<=5; ++i)
    (*tapOut) << setw(8) << "<=" << i << "  ";
  (*tapOut) << setw(8) << ">" << 5 << "  \n";

  // print data line
  const int ROWS = AminoAcid_CODE_SIZE-1;
  int tot = 0;
  for(int i=1; i<=5; ++i)
  {
     tot += getRowNativeCoverage(mtx, chain, start, end, ROWS-i);
     (*tapOut) << setw(9) << tot << "  ";
  }

  tot = 0;
  for(int i=ROWS-6; i>=0; --i)
    tot += getRowNativeCoverage(mtx, chain, start, end, i);
  (*tapOut) << setw(9) << tot << "  \n";

  // delete output file stream, if created
  if(tapOut != &cout)
    delete tapOut;
}



int main(int nArgs, char* argv[])
{
      
  if(nArgs == 1) 
  {
      cerr << "Missing argument specification. Use -h to get help.\n";
      return 1;
  }

  if(getArg("h", nArgs, argv)) 
  {
      sShowHelp();
      return 1;
  };

  if(getArg("-help", nArgs, argv)) 
  {
      sShowSpecial();
      return 1;
  };
  
  // set name of input PDB file
  string inputFile;
  getArg("i", inputFile, nArgs, argv, "!"); 
  if(inputFile == "!") 
  {
    cerr << "Missing input file. Use -h to get help.\n";
    return 1;
  }    

  // set name of knowledge file, default is TOP500
  string know;
  getArg("-know", know, nArgs, argv, "!");
  if(know == "!") 
  {  char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string env;
    env = getenv("VICTOR_ROOT");
    string tmpknow = "data/tor.par";
    string defaultknow = env + tmpknow;
    know = defaultknow;
    if(defaultknow.length() < 3)
	  ERROR("Environment variable VICTOR_ROOT was not found.", exception);
  }

  // select the ensemble of torsional angles used to compute the potential
  bool pp = getArg("-pp", nArgs, argv);      
  bool ppo = getArg("-ppo", nArgs, argv);
  bool ppoc = getArg("-ppoc", nArgs, argv);
  bool ppocc = getArg("-ppocc", nArgs, argv);

  if(((pp) && (ppo)) || ((pp) && (ppoc)) || ((pp) && (ppocc)) ||
      ((ppo) && (ppoc)) || ((ppo) && (ppocc)) || ((ppoc) && (ppocc)))
    ERROR("ERROR in choosing options. You have chosen too many potentials", exception);

  // set bin size for torsional angles, default is 10 degrees
  int ARCSTEP1;
  getArg( "-bin", ARCSTEP1, nArgs, argv, 10);
  if((360%ARCSTEP1) != 0)
    ERROR("The range of angles (phi, psi,) must be divisor of 360.", exception);
  
  // select potential function based on torsional-angle ensemble and bin size
  TorsionPotential* pot = NULL;
  if(pp)
    pot = new PhiPsi(ARCSTEP1, know);
  else if(ppo)
    pot = new PhiPsiOmega(ARCSTEP1, know);
  else if(ppoc)
    pot = new PhiPsiOmegaChi1(ARCSTEP1, know);
  else // default  
    pot = new PhiPsiOmegaChi1Chi2(ARCSTEP1, know); 

  // select chain identifier in the PDB file, default is ' ' 
  string chainID;
  getArg( "-chain", chainID, nArgs, argv, " "); 

  // load protein-chain structure from PDB file
  ifstream inFile(inputFile.c_str());
  if(!inFile)
    ERROR("File " + inputFile + " not found.", exception);
  PdbLoader pl(inFile);

  pl.setChain(chainID.c_str()[0]);

  pl.setNoHAtoms();
  pl.setNoVerbose();
  pl.setPermissive();

  Spacer sp;
  sp.load(pl);
           
  // define the chain fragment [START, END) for which to compute the potential 
  unsigned int START;
  getArg("s", START, nArgs, argv, 9999);
  if(START < 9999)
  {
    START = sp.getIndexFromPdbNumber(START);
    if(START==0)
      ERROR("Impossible to start before the second residue.",exception);
  }
  else
    START = 1;
     
  unsigned int END;
  getArg("e", END, nArgs, argv, 9999);
  if(END < 9999)
    END = sp.getIndexFromPdbNumber(END);
  else
    END = sp.sizeAmino()-1;

  // print program parameters on screen
  bool verbose = getArg("v", nArgs, argv);
  if(verbose)
  {
    cout << "\n#### tap4design ####\n";
    cout << "PDB file: " << inputFile << "\n";
    cout << "chain ID: \"" << chainID << "\"" << "\n";
    cout << "start residue: " << START << "\n";
    cout << "end residue: " << END << "\n";
    cout << "tap function: " << pot->getLabel() << "\n";
    cout << "tap database: " << know << "\n";
    cout << "tap bin size: " << ARCSTEP1 << " degrees.\n\n";
  }

  /////////////////////////////////////////////////

  // build tap profile matrix
  vector<vector<TapRec>*>* tapProfMat = 
    buildFragTapProfileMatrix(sp, START, END, pot);

  // sort position-specific columns by tap score
  const unsigned int fragLen = END - START;
  for(unsigned int i=0; i<fragLen; ++i)
    sortTapRecVec((*tapProfMat)[i]);

  // print matrix and its summary to file
  const char MATRIX_FILE[] = "tap_profile.txt"; 
  printFragTapProfileMatrix3L(tapProfMat, sp, START, END, MATRIX_FILE);  
  printSummaryFTPM(tapProfMat, sp, START, END, MATRIX_FILE);

  // print summary on screen
  printSummaryFTPM(tapProfMat, sp, START, END);

  // delete matrix
  for(unsigned int i=0; i<fragLen; ++i)
    delete (*tapProfMat)[i];
  delete tapProfMat;
}
