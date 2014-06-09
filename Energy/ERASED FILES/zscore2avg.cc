
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <IoTools.h>
#include <StatTools.h>
    
using namespace Biopool;

void sShowHelp()
{
  cout << "Z-score to Average\n"
       << "For FRST Z-score results generated with 'energy2zscore -s'.\n"
       << " Options: \n"
       << "\t-i <filename> \t\t Input file\n"
       << "\t[-v]          \t\t Verbose mode\n"
       << "\t[-l]          \t\t LiveBench/ProQ input data mode\n"
       << endl;
}

int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inputFile;

  getArg( "i", inputFile, nArgs, argv, "!");
  bool verbose = getArg( "v", nArgs, argv);
  bool livebench = getArg( "l", nArgs, argv);

  if (inputFile == "!") 
    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  ifstream inFile(inputFile.c_str());
  if (!inFile)
    ERROR("File not found.", exception);

  unsigned int n_tot = 0;
  unsigned int n_prot = 0;
  unsigned int n_first = 0;
  vector<long double> data, dataCC, dataFE, dataLB1, dataLB10;

  while (inFile)
    {
      string tmp, tmp1, tmp2, tmp3, tmp4;
      int num, rank;
      long double zs, cc = 0.0, fe = 0.0, lb1 = 0.0, lb10 = 0.0;

      if (livebench)
	  inFile >> tmp >> num >> rank >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> zs;
      else
	inFile >> tmp >> num >> zs >> rank >> cc >> fe >> lb1 >> lb10;


      if (!inFile)
	break;
      
      if (verbose)
	cout << tmp << "\n";

      n_tot += num;

      data.push_back(zs);
      dataCC.push_back(cc);
      dataFE.push_back(fe);
      dataLB1.push_back(lb1);
      dataLB10.push_back(lb10);

      if (rank == 1)
	n_first++;
      
      n_prot++;
      skipToNewLine(inFile);
    }

  long double avg = average(data);
  long double sd = standardDeviation(data, avg);
  long double avgCC = average(dataCC);
  long double sdCC = standardDeviation(dataCC, avgCC);
  long double avgFE = average(dataFE);
  long double sdFE = standardDeviation(dataFE, avgFE);
  long double avgLB1 = average(dataLB1);
  long double sdLB1 = standardDeviation(dataLB1, avgLB1);
  long double avgLB10 = average(dataLB10);
  long double sdLB10 = standardDeviation(dataLB10, avgLB10);
  
  cout << "---------------------------------------------------\n";
  
  cout << "N proteins = " << setw(4) << n_prot << "\tN decoys = " 
       << setw(4) << n_tot << "\n"
       << "Average Z-score = " << setw(7) << setprecision(5) << avg 
       << " (+/- " << setw(7) << setprecision(5) << sd << ")\n"
       << "N Rank 1 = " << setw(3) << setprecision(2) << n_first << "\n"
       << "Correlation (Pearson) = " << setw(7) << setprecision(5) << avgCC 
       << " (+/- " << setw(7) << setprecision(3) << sdCC << ")\n"
       << "Frac. Enrich. (10/10) = " << setw(7) << setprecision(5) << avgFE 
       << " (+/- " << setw(7) << setprecision(3) << sdFE << ")\n"
       << "log10 P_best1 =         " << setw(7) << setprecision(5) << avgLB1 
       << " (+/- " << setw(7) << setprecision(3) << sdLB1 << ")\n"
       << "log10 P_best10 =        " << setw(7) << setprecision(5) << avgLB10 
       << " (+/- " << setw(7) << setprecision(3) << sdLB10 << ")\n";
 
  cout << "---------------------------------------------------\n";

  return 0;
}
