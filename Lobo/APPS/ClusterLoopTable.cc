/**
 * @Description This program clusters tables of protein entries for the Loboalgorithm
 */
#include <string>
#include <GetArg.h>
#include <LoopTable.h>
#include <IntCoordConverter.h>

using namespace Biopool;

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))
    cout << "This program clusters tables of protein entries for the Lobo " 
	 << "algorithm.\n"
	 << "Usage: \n" 
	 << " ClusterLoopTable -I<src> -O<file> -C<cutoff> " 
	 << "[-v] [-h] \n" 
	 << "\t -I<file> \t Name of input file.\n"
	 << "\t -O<file> \t Name of the output file.\n"
	 << "\t -C<cutoff> \t Cutoff for clustering. \n"
	 << "\t -v \t\t Verbose.\n"
	 << "\t -h \t\t Display this help.\n"
	 << endl;
  
  bool verbose = getArg( "v", nArgs, argv);
  
  string inputFile, outputFile;
  double cutoff;
  getArg( "I", inputFile, nArgs, argv, "!");
  getArg( "O", outputFile, nArgs, argv, "!");
  getArg( "C", cutoff, nArgs, argv, 0.0);
  
  if ((inputFile == "!") || (outputFile == "!") || (cutoff == 0.0))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  if (verbose)
    cout << "Verify: \n" << "Input from: " << inputFile << "\n" 
	 << "Output to: " << outputFile << "\t Cutoff: " << cutoff << endl;
  
  RamachandranData* ramaData = new RamachandranData;

  LoopTable lt;
  lt.setRama(ramaData);
  lt.read(inputFile);

  lt.cluster(cutoff);

  lt.write(outputFile); 
  
  return 0;
}
