#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <SeqLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <RelSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start\n";
  
  IntCoordConverter icc;
  Spacer sp;

  ifstream refFile("reference.ref");
  ifstream inFile("pdbTest2.seq");

  if (!refFile)
    ERROR("File not found.", exception);
  if (!inFile)
    ERROR("File not found.", exception);

  SeqLoader sl(inFile, refFile);
  sp.load(sl);

  cout << "-------------------------------------------------------\n";
  RelSaver iss(cout);
  sp.save(iss);
  cout << "-------------------------------------------------------\n";
  SeqSaver sss(cout);
  sss.setWriteChi(false);
  sp.save(sss);
  cout << "-------------------------------------------------------\n";

  ofstream pdbFile("pdbRef.ent");

  if (!pdbFile)
    ERROR("File not found.", exception);
  
  PdbSaver pss(pdbFile);
  sp.save(pss);

  cout << "-------------------------------------------------------\n";

  return 0;
}
