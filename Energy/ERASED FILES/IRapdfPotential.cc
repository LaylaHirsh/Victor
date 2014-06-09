/**
 * @Class              IRapdfPotential
 * @Project    Victor
 * @Description 
*    This class implements the *interpolated* all-atom residue specific 
*    probability discriminatory function from Samudrala & Moult (JMB 1998).
*
*/

// Includes:
#include <IRapdfPotential.h>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
string IRapdfPotential::RAPDF_PARAM_FILE = "data/ram.par";

// CONSTRUCTORS/DESTRUCTOR:

IRapdfPotential::IRapdfPotential()
{
     char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
  string inputFile = getenv("VICTOR_ROOT");
  if (inputFile.length() < 3)
    ERROR("Environment variable VICTOR_ROOT was not found.", exception);

  inputFile += RAPDF_PARAM_FILE;

  ifstream in(inputFile.c_str(), ios::in | ios::binary);

  if (!in)
    ERROR("Could not read data file.", exception);

  //load array into heap (ie. faster):
  unsigned int nBins = MAX_BINS * MAX_TYPES * MAX_TYPES;
  double* tmpBin = new double[nBins];
  in.read((char *) (&tmpBin), sizeof(double) * nBins);

  unsigned long compIndex = 0;
  for (unsigned int i = 0; i < MAX_BINS; i++)
    for (unsigned int j = 0; j < MAX_TYPES; j++)
      for (unsigned int k = 0; k < MAX_TYPES; k++)
	{
	  prob[i][j][k] = tmpBin[compIndex];
	  compIndex++;
	}

  in.close();

}

// PREDICATES:
/**
 * @Description Returns the calculated energy
 * @param reference Spacer
 * @return long double
 */
long double IRapdfPotential::calculateEnergy(Spacer& sp)
{
  long double en = 0.0;
  unsigned int size = sp.sizeAmino();

  for (unsigned int i = 0; i < size-1; i++)
      {
	AminoAcid& aa = sp.getAmino(i);
	string aaType = aa.getType();
	for (unsigned int ii = i+1; ii < size; ii++)
	  {
	    AminoAcid& aa2 = sp.getAmino(ii);
	    string aaType2 = aa2.getType();
	    for (unsigned int j = 0; j < aa.size(); j++)
	      for (unsigned int k = 0; k < aa2.size(); k++)
		en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);
	  }
      }
  
  return en;
}
/**
 * @Description Returns the calculated energy
 * @param reference of an AminoAcid and a Spacer
 * @return long double
 */
long double IRapdfPotential::calculateEnergy(AminoAcid& aa, Spacer& sp)
{
  long double en = 0.0;
  string aaType = aa.getType();

  for (unsigned int i = 0; i < sp.sizeAmino(); i++)
    {
      AminoAcid& aa2 = sp.getAmino(i);
      string aaType2 = aa2.getType();
      for (unsigned int j = 0; j < aa.size(); j++)
	for (unsigned int k = 0; k < aa2.size(); k++)
	  en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);
    }
  return en;
}

/**
 * @Description Returns the calculated energy
 * @param reference of two AminoAcids  
 * @return long double
 */
long double IRapdfPotential::calculateEnergy(AminoAcid& aa, AminoAcid& aa2)
{
  long double en = 0.0;
  string aaType = aa.getType();
  string aaType2 = aa2.getType();

  for (unsigned int j = 0; j < aa.size(); j++)
    for (unsigned int k = 0; k < aa2.size(); k++)
      en += calculateEnergy(aa[j], aa2[k], aaType, aaType2);

  return en;
}


unsigned int 
IRapdfPotential::pGetDistanceBinTwo(double distance)
/* No zeros for the distance bin */
{
  if (distance < 11.0)
    {
      if (distance >= 10.0) return 9;
      if (distance >= 9.0) return 8;
      if (distance >= 8.0) return 7;
      if (distance >= 7.0) return 6;
      if (distance >= 6.0) return 5;
      if (distance >= 5.0) return 4;
      if (distance >= 4.0) return 3;
      if (distance >= 3.0) return 2;
      if (distance >= 2.0) return 1;
      return 0; // distance < 3.0
    }
  else
    {
      if (distance >= 18.0) return 17;
      if (distance >= 17.0) return 16;
      if (distance >= 16.0) return 15;
      if (distance >= 15.0) return 14;
      if (distance >= 14.0) return 13;
      if (distance >= 13.0) return 12;
      if (distance >= 12.0) return 11;
      return 10; // distance < 12.0
    }

  ERROR("Unknown distance encountered.", exception);
  return 0;
}
