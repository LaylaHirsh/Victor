/**
 * @Description This program extracts torsion angles from loop regions
 */
#include <string>
#include <GetArg.h>
#include <LoopModel.h>
#include <LoopTableEntry.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <IntCoordConverter.h>
#include <VectorTransformation.h>
#include <limits.h>    
using namespace Biopool;

void sShowHelp(){
  cout << "Loop 2 Torsion\n"
       << "This program extracts torsion angles from loop regions.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input pdb file\n"
       << "\t-c <filename> \t\t chain id\n";
}


void sGetOption(double& param, char* optName, int nArgs, char* argv[], bool verbose = true){
  double tmp = -100000.0;
  getArg(optName, tmp, nArgs, argv, -100000.0);
  if (tmp > -100000.0)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(unsigned int& param, char* optName, int nArgs, char* argv[], bool verbose = true){
  unsigned int tmp = INT_MAX;
  getArg(optName, tmp, nArgs, argv, INT_MAX);
  if (tmp != INT_MAX)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}

void sGetOption(float& param, char* optName, int nArgs, char* argv[], bool verbose = true){
  float tmp = -100000.0;
  getArg(optName, tmp, nArgs, argv, -100000.0);
  if (tmp > -100000.0)    {
      param = tmp;
      if (verbose)
	cout << optName << " = " << tmp << "\n";
    }
}

string sGetSequence(string sequenceFile){
  cout << "seqFile = " << sequenceFile << "\n";
  ifstream seqFile(sequenceFile.c_str());
  if (!seqFile)
    ERROR("Sequence file not found.", exception);

  string tmp;
  seqFile >> tmp;
  return tmp;
}

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile,chainID;
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, " ");
  if ((inputFile == "!"))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  Spacer *sp;

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile);
  pl.setNoHAtoms();
  vector<char> allCh; 
   allCh = pl.getAllChains(); 
      for (unsigned int i = 0; i < allCh.size(); i++)
           cout << "\t," << allCh[i] << ",";
      cout << "\n";

    /*check on validity of chain: 
    if user select a chain then check validity
     else select first valid one by default*/
      if (chainID != " ")  {
        bool validChain=false;
        for (unsigned int i = 0; i < allCh.size(); i++ ) {
          if (allCh[i]==chainID[0])  {
            pl.setChain(chainID[0]);
            cout << "Loading chain " << chainID << "\n";
            validChain=true;
            break;
          }
        }
        if (!validChain) {
           cout << "Chain " << chainID << " is not available\n";
          return -1;
        }

      }
    
  Protein prot;
  prot.load(pl);
  sp=prot.getSpacer(chainID[0]);
  

  sp->setStateFromTorsionAngles();
 
  for (unsigned int i = 1; i < sp->sizeAmino()-1; i++)    {
      if (sp->getAmino(i).getState() == COIL)
	cout << "   " << setw(5) << setprecision(3) << sp->getAmino(i).getPhi()
	     << "   " << setw(5) << setprecision(3) << sp->getAmino(i).getPsi()
	     << "    1.0\n";
    }


  return 0;
}
