/**
 * @Description This program clusters the data contained in a Ramachandran distribution file
 */
#include <string>
#include <GetArg.h>
#include <RamachandranData.h>

using namespace Biopool;

int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))
    cout << "This program clusters the data contained in a Ramachandran"
	 << " distribution file.\n"
	 << "Usage: \n" 
	 << " LoopTableTest -A<src1> -B<src2> -O<file> -R<rama> -S<size> " 
	 << "[-v] [-h] \n" 
	 << "\t -i<file> \t Name of input file.\n"
	 << "\t -o<file> \t Name of the output file.\n"
	 << "\t -c<size> \t Cutoff value.\n"
	 << "\t -h \t\t Display this help.\n"
	 << endl;
  
  string inputFile, outputFile;
  double cutoff;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "c", cutoff, nArgs, argv, 0.0);
  
  if ((inputFile == "!") || (outputFile == "!") || (cutoff == 0.0))     {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }
  
  RamachandranData* ramaData = new RamachandranData;
  ifstream ramaIn(inputFile.c_str());
  if (!ramaIn)
    ERROR("Rama filename invalid.", exception);
  ramaData->load(ramaIn);
  
  ramaData->cluster(cutoff);

  ofstream ramaOut(outputFile.c_str());
  if (!ramaOut)
    ERROR("Rama filename invalid.", exception);
  ramaData->save(ramaOut);

  return 0;
}
