#include <Spacer.h>
#include <Group.h>
#include <SideChain.h>
#include <vector3.h>
#include <SeqLoader.h>
#include <PdbSaver.h>
#include <XyzSaver.h>
#include <SeqSaver.h>
#include <IoTools.h>
#include <IntCoordConverter.h>
#include <SeqConstructor.h>
#include <IntSaver.h>

using namespace Biopool;

int main()
{ 
  cout << "Start" << endl;

  IntCoordConverter icc;

  Spacer sp;

  ifstream refFile("reference.ref");

  if (!refFile)
    ERROR("File not found.", exception);

  SeqConstructor sc(refFile);

  cout << "-------------------------------------------------------\n";

  sp = sc.makeSpacer("GLY", 12);

  cout << "-------------------------------------------------------\n";
  cout << "Writing XYZ file...\n";

  ofstream outFile("builder.xyz");

  if (!outFile)
    ERROR("File not found.", exception);

  XyzSaver xs(outFile);
  sp.save(xs);
  cout << "-------------------------------------------------------\n";
  cout << "Writing PDB file...\n";

  ofstream outFile2("builder.pdb");

  if (!outFile2)
    ERROR("File not found.", exception);

  PdbSaver ps(outFile2);
  sp.save(ps);
  cout << "-------------------------------------------------------\n";
  SeqSaver is(cout);
  sp.save(is);
  cout << "-------------------------------------------------------\n";


  return 0;
}
