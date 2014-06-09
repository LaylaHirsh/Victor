/**
 * @Description This program generates GNU compatible plots of protein entries  for the Lobo algorithm.
 */
#include <string>
#include <LoopTable.h>
#include <stdlib.h>
#include <GetArg.h>
#include <VectorTransformation.h>

using namespace Biopool;

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))    {
      cout << "Loop Table Plot\n"
	   << "This program generates GNU compatible plots of protein entries" 
	   << " for the Lobo algorithm.\n"
	   << " Options: \n"
	   << "\t-i <filename> \t\t Input filename\n"
	   << "\t-o <filename> \t\t Output filename\n"
	   << "\t-s <s m l> \t\t Size: small, medium, large\n"
	   << "\t[-w <0 1 2>] \t\t Which vector: EP, EN, or ED"
	   << " (default= 0)\n"
	   << "\t[-d <xy/xz/yz: 0 1 2>] \t Which direction: xy, xz or yz" 
	   << " (default= 0)\n"
	   << "\t[--noRot] \t\t Do not rotate table into xy plane\n"
	   << "\t[--show] \t\t Show composition of distance bins.\n";
      return 1;
    };

  string inputFile, outputFile, size;
  unsigned int selVec, selDir;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "s", size, nArgs, argv, "!");
  getArg( "w", selVec, nArgs, argv, 0);
  getArg( "d", selDir, nArgs, argv, 0);
  bool noRot = getArg( "-noRot", nArgs, argv);
  bool show = getArg( "-show", nArgs, argv);
 
 if ((inputFile == "!") || (outputFile == "!"))    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  cout << "Verify: \n" << "Input from: " << inputFile << "\n" 
       << "Output to: " << outputFile << "\t Size = " << size << endl;
  
  LoopTable P;
  VectorTransformation vt;
  
  P.read(inputFile);

  if (show)
    P.showDistribution();

  if (!noRot)
    for (unsigned int j = 0; j < P.size(); j++)
      P[j].rotateIntoXYPlane(vt);

  if (size == "s")
    P.writeASCII(outputFile, 1000, selVec, selDir);
  else if (size == "m")
    P.writeASCII(outputFile, 10000, selVec, selDir);
  else if (size == "l")
    P.writeASCII(outputFile, 100000, selVec, selDir);
  else
    cout << "Error selecting table size." << endl;

  return 0;
}
