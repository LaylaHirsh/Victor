/**  
@Description */
#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <XyzLoader.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Biopool;

int main(int nArgs, char* argv[])
{ 
  if (nArgs != 2)
    {
      cout << "Pdb Secondary $Revision: 1.2 $ -- output:\n" 
	   << "1) aminoacid sequence (one letter code);\n"
	   << "2) the type of aminoacid " 
	   << "\t N = negative charge, \t P = positive charge,"
	   << " \t h = hydrophilic, \t + = hydrophobic, "
	   << "\t , = neutral (hydrophobic);\n"
	   << "3) the states of all aminoacids in the protein " 
	   << "\t H = helix, \t E = extended (strand, sheet), \t . = other.\n"
	   << "\n";
      return 1;
    };

  Spacer sp;
  ifstream inFile(argv[1]);
  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader il(inFile);
  sp.load(il); 

  for (unsigned int i = 1; i <= sp.sizeAmino(); i++)
    {
      if (i % 10 == 0)
	cout << setw(1) << ((i%100)/10);
      else 
	cout << " ";
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    {
      cout << sp.getAmino(i).getType1L();
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";
  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    {
      switch (sp.getAmino(i).getCode())
	{
	case ASP:
	case GLU:
	  cout << "P"; break;
	case LYS:
	case ARG:
	  cout << "N"; break;
	case ASN:
	case GLN:
	case SER:
	case THR:
	case HIS:
	  cout << "h"; break;
	case VAL:
	case LEU:
	case ILE:
	  cout << "+"; break;
	default:
	  cout << ",";
      };
      if ((i+1) % 60 == 0)
	cout << "\n";
    }
  cout << "\n";
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
