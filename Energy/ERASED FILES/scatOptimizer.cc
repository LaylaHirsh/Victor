/**  
@Description */
#include <EnergyAnalyzer.h>
#include <GetArg.h>
#include <iostream>
#include <IoTools.h>
//#include <StatTools.h>

using namespace Biopool;

void sShowHelp()
{
  cout << "Scat Edit\n"
       << "This program allows to edit scatter plots from MQAP programs.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t[-o <filename>] \t Output file (default = scat.plot)\n"
       << "\t-x <column> \t\t Select column for x axis in output\n"
       << "\t-y <column> \t\t Select column for y axis in output\n"
       << "\t\t NB: use code 99 to select optimized score I in -x & -y \n"
       << "\t [--noOpt1] \t do not optimize with scheme I\n"
       << "\t [-m] \t multiply GDT_TS score by 100 I\n"; 
}

void sLine()
{
  cout << "-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n";
}

double sCorrelation(vector<double>& zs1,vector<double>& zs2 )
{
    double tmp = 0.0;
    for (unsigned int j = 0; j < zs1.size(); j++)
	tmp += zs1[j] * zs2[j];

    tmp /= zs1.size();
    return tmp;
}


int main(int nArgs, char* argv[])
{ 
  cout.setf(ios::fixed, ios::floatfield);

  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile;

  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "scat.plot");

  unsigned int colX,colY;
  getArg( "x", colX, nArgs, argv, 0);
  getArg( "y", colY, nArgs, argv, 0);

  if (inputFile == "!") 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  bool noOpt1 = getArg( "-noOpt1", nArgs, argv);
  bool mult = getArg( "m", nArgs, argv);

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  EnergyAnalyzer ra;

  ra.load(inFile, 0, mult);
  ra.calcStatistics();

  // correlations:
  ra.printAllCorrelations();

  vector<double> bestRaw;
  // optimization:
  if (!noOpt1)
    bestRaw = ra.linearOptimization(0.5, 3.5, 0.125, 0.0, 3000.0, 50.0,
		       -200.0, 200.0, 25.0, 0.0, 3000.0, 50.0);
//   else
//     bestRaw = ra.nonLinearOptimization(100, 100);

  // write plot file (optional):
  if (colX * colY != 0)
  {
    cout << "Writing scatter plot to " << outputFile << "\n";
      ofstream out(outputFile.c_str());
      if (!out)
	  ERROR("File not found.", exception);

      for (unsigned int i = 0; i < ra.col[0].size(); i++)
      {
	if (colX == 99)
	  out << setw(5) << bestRaw[i] << "\t";
	else
	  out << setw(5) << ra.col[colX-1][i] << "\t";

	if (colY == 99)
	  out << setw(5) << bestRaw[i] << "\n";
	else
	  out << setw(5) << ra.col[colY-1][i] << "\n";
      }
      sLine();
  }

    return 0;
}
