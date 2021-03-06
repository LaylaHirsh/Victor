/**
 @Description FRST - Function of Rapdf, Solvation and Torsion potentials Structure refinement by torsion angle minimization."
 */ 
#include <string>
#include <GetArg.h>
#include <PdbLoader.h>
#include <PdbSaver.h>
#include <SolvationPotential.h>
#include <RapdfPotential.h>
#include <TorsionPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmega.h>
#include <PhiPsiPreAngle.h>
#include <PhiPsiOmegaChi1Chi2PreAngle.h>

using namespace Biopool;

const long double W_RAPDF = 2.5;
const long double W_SOLV = 500.0;
const long double W_HYDB = -50.0;
const long double W_TORS = 350.0;

void sShowHelp(){
  cout << "FRST - Function of Rapdf, Solvation and Torsion potentials\n\n"
       << "Structure refinement by torsion angle minimization."
       << "Developed by Silvio Tosatto <silvio.tosatto@unipd.it> \n"
       << "   Options: \n"
       << "\t-i <filename> \t\t Input file for PDB template\n"
       << "\t-o <filename> \t\t Output file for refined PDB structure\n"
       << "\t[-c <id>]  \t ID of chain to load from PDB file\n"
       << "\t[-r <integer>] \t\t Maximum number of refinement iterations (default = 1)\n"
       << "\t[-m <integer>] \t\t Maximum number of refinement steps (default = 50)\n"
       << "\t[-s <integer>] \t\t Maximum local torsion angle step (default = 2)\n"
       << "\t[-t] \t\t\t New (Mk2) torsion angle potential \n"
       << "\t[-v] \t\t\t Verbose mode\n"
       
       << "\n";
}

double sClashScore(Spacer& sp, unsigned int pos){
  double clash = 0.0;

  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    if (i != pos)      {
	if (sp.getAmino(i)[CA].distance(sp.getAmino(pos)[CA]) < 3.2) {
	    clash += 1;
	    cout << "I";
	  }
      }
  return clash;
}

double sHydrogen(Spacer& sp){
  int count = 0;
  
  for (int i = 0; i < sp.sizeAmino(); i++)
    for (int j = 0; j < sp.sizeAmino(); j++)
      if (abs(i-j) > 1)	{
	  double dist = sp.getAmino(i)[N].distance(sp.getAmino(j)[O]);
	  double dist2 = sp.getAmino(i)[N].distance(sp.getAmino(j)[C]);
	  double dist3 = sp.getAmino(i)[CA].distance(sp.getAmino(j)[O]);
	  
	  if ((dist <= 4.0) && (dist >= 2.0) && (dist < dist2) && (dist < dist3))   {
	      count++;
	      break;
	    }
	}
  return count;
}


int main(int nArgs, char* argv[]){ 
  if (getArg( "h", nArgs, argv))  {
      sShowHelp();
      return 1;
    };

  string inputFile, outputFile, chainID;
  unsigned int MAX_REF, maxIter;
  double step;
  bool newTorsion = getArg( "t", nArgs, argv);
  bool verbose = getArg( "v", nArgs, argv);
  getArg( "i", inputFile, nArgs, argv, "!");
  getArg( "o", outputFile, nArgs, argv, "!");
  getArg( "c", chainID, nArgs, argv, " ");
  getArg( "m", MAX_REF, nArgs, argv, 50);
  getArg( "r", maxIter, nArgs, argv, 1);
  getArg( "s", step, nArgs, argv, 2);
  
  if ((inputFile == "!") || (outputFile == "!"))    {
      cout << "Missing file specification. Aborting. (-h for help)" << endl;
      return -1;
    }

  SolvationPotential solv;
  RapdfPotential rapdf;

  TorsionPotential* tors = NULL;
  if (!newTorsion)
    tors =  new PhiPsi(5);//default ARCSTEP = 5
  else
    tors =  new PhiPsiOmegaChi1Chi2PreAngle(20);//ARCSTEP 20, ARCSTEP2 40 
   
  Spacer *sp;
  ifstream inFile2(inputFile.c_str());
  if (!inFile2)
    ERROR("File not found.", exception);
  PdbLoader pl(inFile2);
  pl.setChain(chainID.c_str()[0]);
  pl.setNoHAtoms();
  pl.setNoVerbose();
  
  pl.setPermissive();
  Protein prot;
  prot.load(pl);
  sp=prot.getSpacer(chainID.c_str()[0]);
  if (!pl.isValid())
      ERROR("Invalid PDB file found.\n", exception);

  cout.setf(ios::fixed, ios::floatfield);
   
  for (unsigned int g = 0; g < maxIter; g++)
       for (unsigned int i = 0; i < sp->sizeAmino()-1; i++)    {
      unsigned int pos = i;

      cout << "pos= " << setw(4) << pos << "   ";

      double phi = sp->getAmino(pos).getPhi();
      double psi = sp->getAmino(pos).getPsi();
      double omega = sp->getAmino(pos).getOmega();

      cout << "( " << setw(7) << setprecision(2) << phi << " / " 
	   << setw(7) << setprecision(2) << psi << " / " 
	   << setw(7) << setprecision(2) << omega << ") \t";

      double torsEn = dynamic_cast<PhiPsi*>(tors)->calculateEnergy(sp->getAmino(pos)) 
	+ sClashScore(*sp, pos);

      cout << "en = " << setw(7) << setprecision(2) << torsEn << "\t";

      // search for gradient:
      double phiOpt = phi;
      double psiOpt = psi;
      double omegaOpt = omega;

      double l = 0.0;
      for (double j = -step; j <= step; j += step/4)
	for (double k = -step; k <= step; k += step/4)	  {
	    sp->getAmino(pos).setPhi(phi + j);
	    sp->getAmino(pos).setPsi(psi + k);
	    sp->getAmino(pos).setOmega(omega + l);
	    sp->getAmino(pos).sync();

	    double tmp = dynamic_cast<PhiPsi*>(tors)->calculateEnergy(sp->getAmino(pos)) 
	      + sClashScore(*sp, pos);
	    if (tmp < torsEn)	      {
		torsEn = tmp;
		phiOpt = phi + j;
		psiOpt = psi + k;
		omegaOpt = omega + l;
	    }
	  }

      sp->getAmino(pos).setPhi(phiOpt);
      sp->getAmino(pos).setPsi(psiOpt);
      sp->getAmino(pos).setOmega(omegaOpt);
      sp->getAmino(pos).sync();

      cout << "--> ( " << setw(7) << setprecision(2) << phiOpt << " / " 
	   << setw(7) << setprecision(2) << psiOpt << " / " 
	   << setw(7) << setprecision(2) << omegaOpt << ") \t" 
	   << "en = " << setw(7) << setprecision(2) << torsEn;
	    
      cout << "\n";
    }

  ofstream outFile2(outputFile.c_str());
  if (!outFile2)
    ERROR("File not found.", exception);
  PdbSaver ps(outFile2);
  sp->save(ps);

  delete tors;
}
