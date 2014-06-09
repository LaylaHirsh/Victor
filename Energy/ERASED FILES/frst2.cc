/**
 @Description FRST2 - Function of Rapdf, Solvation and Torsion potentials, Mk2
 * This program calculates several pseudo-energies and features 
 * to evaluate the quality of a given protein structural model with
 * a machine learning approach (epsilon-SVR).  
 */ 
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <EnergyFeatures.h>
#include "libsvm_wrapper.h"

using namespace Biopool;

void sShowHelp()
{
  cout <<"FRST2 - Function of Rapdf, Solvation and Torsion potentials, Mk2\n\n"
       << "This program calculates several pseudo-energies and features "
       << "to evaluate the quality of a given protein structural model with"
       << "a machine learning approach (epsilon-SVR). Contact the author for "
       << "further details.\n"
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n"
       << "   Options: \n"
       << "\t-I <filelist> \t\t Input *filelist* file for PDB templates\n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t-m <filename> \t\t SVM model file (libsvm format)\n"
       << "\t-r <filename> \t\t SVM features range file (libsvm format)\n"
       << "\t[--chain <id>]  \t ID of chain to load from PDB file\n"
       << "\t[--features]  \t Shows the list of features\n"
       << "\t[--quiet]  \t Skip target name (useful for benchmarking)\n"
       << "\t[--verbose]  \t Print individual features\n"
       << "\n";
}


int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile, inputFilelist, chainID, modelfname, rangefname;
  bool showFeatures = getArg( "-features", nArgs, argv);
  bool quiet = getArg( "-quiet", nArgs, argv);
  bool verbose = getArg( "-verbose", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "I", inputFilelist, nArgs, argv, "!");
  getArg( "m", modelfname, nArgs, argv, " ");
  getArg( "r", rangefname, nArgs, argv, " ");
  getArg( "-chain", chainID, nArgs, argv, " ");

  if ((modelfname == " ") || (rangefname == " "))
    { char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
      string victorRoot = getenv("VICTOR_ROOT");
      if (victorRoot.length() < 3)
	ERROR("Environment variable VICTOR_ROOT was not found.", exception);
      
      if (modelfname == " ")
	modelfname = victorRoot + "data/frst2.model";
      if (rangefname == " ")
	rangefname = victorRoot + "data/frst2.range";
    }

  if (showFeatures)
    {
      EnergyFeatures::showFeatures();
      return 1;
    }

  if ((inputFile == "!") && (inputFilelist == "!"))
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  if ((inputFile != "!") && (inputFilelist != "!"))
    {
       cout << "Please choose between filelist and file mode. Aborting. "
	    << "(-h for help)" << endl;
      return -2;     
    }
    
  EnergyFeatures feat;
  LibSVMWrapper predictor(modelfname, rangefname);
  
  ifstream inFile(inputFilelist.c_str());
  if ((!inFile) && (inputFilelist != "!"))
    ERROR("File not found.", exception);

  while ((inFile) || (inputFile != "!"))
    {
      if (inputFilelist != "!")
	{
	  inFile >> inputFile;
	  if (!inFile)
	    break;
	}

      Spacer sp;
      ifstream inFile2(inputFile.c_str());
      if (!inFile2)
	ERROR("File not found.", exception);
      PdbLoader pl(inFile2);
      pl.setChain(chainID.c_str()[0]);
      pl.setNoHAtoms();
      pl.setNoVerbose();

      pl.setPermissive();
      sp.load(pl);

      if (!pl.isValid())
	{
	  cerr << "Warning: Invalid PDB file found.\n";
	  inputFile = "!";
	  continue;
	}

      vector<double> tmpF = feat.calculateFeatures(sp);

      cout.setf(ios::fixed, ios::floatfield);

      if (!quiet)
	cout << inputFile << "\t";
      cout << predictor(tmpF);

      if (verbose)
	{
	  for (unsigned int i = 0; i < tmpF.size(); i++)
	    cout  << "\t" << setw(9) << setprecision(4) << tmpF[i];
	}
      cout << endl;

      // reset variable to trigger break condition:
      inputFile = "!";
   }
  
}
