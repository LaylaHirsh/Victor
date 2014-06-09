/**
*/
#include <Spacer.h>
#include <PdbLoader.h>
#include <SeqSaver.h>
#include <IoTools.h>

using namespace Biopool;

int main(int nArgs, char* argv[])
{ 
  if (nArgs != 3)
    {
      cout << "Pdb 2 Seq $Revision: 1.2 $ -- converts a PDB file into SEQ " 
	   << "(torsion angles) protein structure backbone torsion angles" 
	   << endl;
      cout << "  Usage: \t\t Pdb2Seq <input_filename> <output_filename> \n";
      return 1;
    };

  Spacer sp;
  ifstream inFile(argv[1]);
  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader il(inFile);
  il.setNoHAtoms();
  sp.load(il); 

  ofstream outFile(argv[2]);
  if (!outFile)
    ERROR("Couldn't write file.", exception);
  
  SeqSaver ss(outFile);
  sp.save(ss);

  return 0;
}
