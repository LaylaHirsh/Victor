#include <AminoAcid.h>
#include <Group.h>
#include <vector3.h>
#include <XyzLoader.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start" << endl;

  AminoAcid aa;
  ifstream inFile("amino.xyz");
  if (!inFile)
    ERROR("File not found.", exception);

  XyzLoader il(inFile);
  aa.load(il);

  cout << "-------------------------------------------------------\n";
  XyzSaver is(cout, 0);
  aa.save(is);

  cout << "-------------------------------------------------------\n";
  cout << "SideChain->backbone = " << aa.getSideChain().getBackboneRef()->getType() << "\n";

  cout << "-------------------------------------------------------\n";
  SeqSaver iss(cout);
  aa.save(iss);
  cout << "-------------------------------------------------------\n";
  cout << "superior: " << "\t aa[2]= " << aa[2].getSuperior().getType()
       << " (" << aa[2].getSuperior().getCode() << ") \t aa[8]= " 
       << aa[8].getSuperior().getType() 
       << " (" << aa[8].getSuperior().getCode() << ")" 
       << endl;

  cout << "-------------------------------------------------------\n";
  IntSaver intS(cout);
  aa.save(intS);
  cout << "-------------------------------------------------------\n";
  return 0;
}
