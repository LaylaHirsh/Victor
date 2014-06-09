/**  
@Description */

#include <string>
#include <GetArg.h>
#include <Spacer.h>
#include <PdbLoader.h>
#include <IoTools.h>

using namespace Biopool;
void sShowHelp()
{
  cout << "Pdb 2 Secondary Structure converter\n" 
       << "\t H = helix, \t E = extended (strand, sheet), \t . = other.\n"
       << "   Options: \n"
       << "\t-i <filename> \t\t Input file for PDB structure\n"
       << "\n";
}

int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile;
  getArg( "i", inputFile, nArgs, argv, "!");

  if (inputFile == "!")
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  Spacer sp;
  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader il(inFile);
  il.setNoHAtoms();
  sp.load(il); 

  cout << ">" << inputFile << "\n";
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    {
      switch (sp.getAmino(i).getState()) 
	{
	case HELIX:
	  cout << "H"; break;
	case STRAND:
	  cout << "E"; break;
	default:
	  cout << ".";
	};
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";
  
  return 0;
}
