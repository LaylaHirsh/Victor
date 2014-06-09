/**  
@Description */

#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <EffectiveSolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2.h>
#include <AminoAcidCode.h>

using namespace Biopool;

void sShowHelp()
{
  cout << "sum4pdb\n\n"
       << "Summarizes per-chain data into per-protein data. To use with WHAT_CHECK.\n"
       << "Options:\n"
       << "\t -i <filename> \t Input file\n"
       << "\t -o <filename> \t Output file\n"
       << endl;
}
int main(int nArgs, char* argv[])
{ 
  if (getArg( "h", nArgs, argv))
    {
      sShowHelp();
      return 1;
    };

  string inFile, outFile;
  getArg( "i", inFile, nArgs, argv, "!");
  getArg( "o", outFile, nArgs, argv, "!");

  ifstream in(inFile.c_str());
  if (!in)
    ERROR("File not found.", exception);
  
  ofstream out(outFile.c_str());
  if (!out)
    ERROR("Could not open file for writing.", exception);

  string code, qual, tors, phipsi, chi, back, rot, chi12, pack1, pack2;
  string prevPdb = "";

  double qu = 0.0; 
  double tor = 0.0;
  double phi = 0.0;
  double ch = 0.0;
  double ba = 0.0;
  double ro = 0.0;
  double c1 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;
  unsigned int countQu = 0;
  unsigned int countTor = 0;
  unsigned int countPhi = 0;
  unsigned int countCh = 0;
  unsigned int countBa = 0;
  unsigned int countRo = 0;
  unsigned int countC1 = 0;
  unsigned int countP1 = 0;
  unsigned int countP2 = 0;

  while (in)
    {
      in >> code >> qual >> tors >> phipsi >> chi >> back >> rot >> chi12 
	 >> pack1 >> pack2;

      if (qual == "na")
	continue;

      string pdb = code.substr(0,4);

      if ((pdb != prevPdb) && (prevPdb != "")) 
	{
	  out << prevPdb 
	      << " \t " << setw(6) << setprecision(5) << qu/countQu 
	      << " \t " << setw(6) << setprecision(5) << tor/countTor 
	      << " \t " << setw(6) << setprecision(5) << phi/countPhi 
	      << " \t " << setw(6) << setprecision(5) << ch/countCh 
	      << " \t " << setw(6) << setprecision(5) << ba/countBa 
	      << " \t " << setw(6) << setprecision(5) << ro/countRo 
	      << " \t " << setw(6) << setprecision(5) << c1/countC1 
	      << " \t " << setw(6) << setprecision(5) << p1/countP1 
	      << " \t " << setw(6) << setprecision(5) << p2/countP2 
	      << "\n";
	  qu = 0.0; 
	  tor = 0.0;
	  phi = 0.0;
	  ch = 0.0;
	  ba = 0.0;
	  ro = 0.0;
	  c1 = 0.0;
	  p1 = 0.0;
	  p2 = 0.0;
	  countQu = 0;
	  countTor = 0;
	  countPhi = 0;
	  countCh = 0;
	  countBa = 0;
	  countRo = 0;
	  countC1 = 0;
	  countP1 = 0;
	  countP2 = 0;
	}

      prevPdb = pdb; 

      if (qual != "na")
	{
	  qu += stod(qual);
	  countQu++;
	}

      if (tors != "na")
	{
	  tor += stod(tors);
	  countTor++;
	}

      if (phipsi != "na")
	{
	  phi += stod(phipsi);
	  countPhi++;
	}

      if (chi != "na")
	{
	  ch += stod(chi);
	  countCh++;
	}

      if (back  != "na")
	{
	  ba += stod(back);
	  countBa++;
	}

      if (rot != "na")
	{
	  ro += stod(rot);
	  countRo++;
	}

      if (chi12 != "na")
	{
	  c1 += stod(chi12);
	  countC1++;
	}

      if (pack1 != "na")
	{
	  p1 += stod(pack1);
	  countP1++;
	}

      if (pack2 != "na")
	{
	  p2 += stod(pack2);
	  countP2++;
	}

    }

}
