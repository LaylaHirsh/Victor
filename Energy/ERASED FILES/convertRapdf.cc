/**  
@Description */
#include <string>
#include <GetArg.h>
#include <RapdfPotential.h>

using namespace Biopool;

void sShowHelp()
{
  cout << "Convert RAPDF Data File to Binary\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file (RAPDF ASCII)\n"
       << "\t-o <filename> \t\t Output file (RAPDF BINARY)\n"
       << "\n";
}


int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile;
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "test.out");

  if ( (inputFile == "!") || (outputFile == "!") )
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

//    RapdfPotential rapdf;

//    rapdf.readASCII(inputFile);
//    cout << "loaded." << endl;

//    rapdf.write(outputFile);
//    cout << "saved." << endl;


//  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";

}
