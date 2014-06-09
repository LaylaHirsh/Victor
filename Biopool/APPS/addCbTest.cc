/**  
@Description      A simple program to test class Protein's features.*/
#include <GetArg.h>
#include <Protein.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>
#include <iostream>
using namespace Biopool;

void sShowHelp()
{
  cout << "addCb test\n"
       << "Options: \n"
       << "\t-i <filename> \t Input PDB file\n"
       << "\t-o <filename> \t Output to file\n"
       << "\n";
}

void sAddLine()
{
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

int main(int argc, char* argv[])
{ 
  if (getArg( "h", argc, argv))
  {
      sShowHelp();
      return 1;
  }
  string inputFile,outputFile,chainID;
  
  getArg( "i", inputFile, argc, argv, "!");
  getArg( "o", outputFile, argc, argv, "!");
  if (inputFile == "!")
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);     
  
  PdbLoader pl(inFile);
  
  pl.setVerbose();
  pl.setNoSecondary();
  pl.setNoHAtoms();                     //noHAtoms=true
  pl.setNoHetAtoms();               //noHetAtoms=true
  
  Spacer sp;
  sp.load(pl);
  cout<<"Spacer Loaded\n"
      <<"...\n";
  
  //section to check what loaded
  cout<<"Spacer size: "<<sp.sizeAmino()<<endl;
  cout<<"inizio patching Cb...\n";
  for(unsigned int i=0; i<sp.sizeAmino(); i++)
  {
	  AminoAcid& am=sp.getAmino(i);
	  cout<<"Tipo: "<<am.getType1L()<<endl;
	  if (am.getSideChain().isMember(CB))
	  {
		  cout<<"Amino "<<i<<" has Cb carbon"<<endl;
	  }
	  else if (am.getType()!= "GLY")
	  {
		  cout<<"Adding Cb carbon to Ammino"<<i<<endl;
		  am.patchBetaPosition();
	  }
	  else
	  {
		  cout<<"Ammino is a Glycina. No Cb in this Case\n";
	  }
  }
  if (outputFile == "!")
    {
      cout << "Missing output file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  ofstream os(outputFile.c_str());
  if (!os)
    ERROR("Could not open file for writing.", exception);
  PdbSaver ps(os);
  ps.setWriteAtomOnly();
  sp.save(ps);
  cout<<"...\n"
      <<"Spacer saved"<<endl;
  
  sAddLine();
  cout<<"EXIT addCbTest" << endl;
  return 0;
}

