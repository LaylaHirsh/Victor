/**

*/
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <IoTools.h>

using namespace Biopool;

int main(int nArgs, char* argv[])
{ 
  if (nArgs != 3)
    {
      cout << "Pdb Editor $Revision: 1.2 $ -- allows sequential manipulation of "
	   << "protein structure backbone torsion angles" << endl;
      cout << "  Usage: \t\t PdbEditor <input_filename> <output_filename> \n";
      return 1;
    };

  Spacer sp;
  ifstream inFile(argv[1]);
  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader il(inFile); 
  sp.load(il); 

  cout << "Editing " << argv[1] << " output goes to " << argv[2] << "\n";

  int aaid = -1;
  do {
    cout << "Aminoacid# or -1: ";
    cin >> aaid;
    if (aaid <= -1)
      {
	cout << "Bye.\n";
	return 0;
      };
    
    if (aaid >= (int)sp.sizeAmino())
      {
      cout << "\t Invalid aa#!\n";
      }
    else
      {
	double newVal = 999;
	cout << "\t " << aaid << "  " << sp.getAmino(aaid).getType() << "\n";
	cout << "Phi=   " << sp.getAmino(aaid).getPhi() 
	     << "\t new phi or 999:   ";
	cin >> newVal;
	if (newVal != 999)
	  sp.getAmino(aaid).setPhi(newVal);
	cout << "Psi=   " << sp.getAmino(aaid).getPsi() 
	     << "\t new psi or 999:   ";
	cin >> newVal;
	if (newVal != 999)
	  sp.getAmino(aaid).setPsi(newVal);
	cout << "Omega= " << sp.getAmino(aaid).getOmega() 
	     << "\t new omega or 999: ";
	cin >> newVal;
	if (newVal != 999)
	  sp.getAmino(aaid).setOmega(newVal);

	ofstream outFile2(argv[2]);
	
	if (!outFile2)
	  ERROR("Couldn't write file.", exception);
	
	PdbSaver pss2(outFile2);
  
	sp.save(pss2);
      };
  }
  while (aaid != -1);
  
  return 0;
}
