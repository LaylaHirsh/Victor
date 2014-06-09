
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
      cout << "Pdb Corrector $Revision: 1.2 $ -- adds missing oxygen atoms to " 
	   << "protein structure backbones" << endl;
      cout << "  Usage: \t\t PdbCorrector <filename> \n";
      return 1;
    };

  Spacer sp;
  ifstream inFile(argv[1]);
  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader il(inFile);
  sp.load(il); 

  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    sp.getAmino(i).addMissingO();

  ofstream outFile2(argv[1]);

  if (!outFile2)
    ERROR("Couldn't write file.", exception);

  PdbSaver pss2(outFile2);
  sp.save(pss2);
  
  return 0;
}
