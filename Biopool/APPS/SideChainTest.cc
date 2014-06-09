/**  
@Description */
#include <Spacer.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SideChainConstructor.h>
#include <GetArg.h>

using namespace Biopool;

int main(int nArgs, char* argv[])
{ 
  cout << "Start" << endl;

  string inputFile, outputFile;
  getArg( "i", inputFile, nArgs, argv, "test.pdb");
  getArg( "o", outputFile, nArgs, argv, "sctest.pdb");

  IntCoordConverter icc;

  Spacer sp;

  ifstream inFile(inputFile.c_str());

  if (!inFile)
    ERROR("File not found.", exception);

  PdbLoader pl(inFile);
  sp.load(pl);

  ifstream refFile("allaminoacids.pdb");

  if (!refFile)
    ERROR("File not found.", exception);

  SideChainConstructor scc(refFile);
  vector<SideChain> vsc;
  
  cout << "-------------------------------------------------------\n";

  for (unsigned int i = 0; i < sp.size(); i++)
    sp.getAmino(i).setSideChain(scc.makeSideChain(sp.getAmino(i).getType()));
  
  cout << "-------------------------------------------------------\n";
  ofstream outFile(outputFile.c_str());

  if (!outFile)
    ERROR("File not found.", exception);

  PdbSaver ps(outFile);
  sp.save(ps);
  cout << "-------------------------------------------------------\n";


  return 0;
}
