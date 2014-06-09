/**
 * @Class               MaxPropensitiesYPreAngle
 * @Project        Victor
*/
// Includes:
#include <MaxPropensitiesYPreAngle.h>

using namespace Biopool;

// Global constants, typedefs, etc. (to avoid):
 string MaxPropensitiesYPreAngle::TOR_PARAM_FILE = "data/tor.par";
 int MaxPropensitiesYPreAngle::ARC_STEP = 10;
 int MaxPropensitiesYPreAngle::SIZE_OF_TABLE =360/ARC_STEP;

// CONSTRUCTORS/DESTRUCTOR:

MaxPropensitiesYPreAngle::MaxPropensitiesYPreAngle() :
   PARAM_FILE ("/data/EnergiaMinAA.txt")
{
  pConstructData();
}

void MaxPropensitiesYPreAngle::pConstructData()
{
   char *victor=getenv("VICTOR_ROOT");
	if (victor  == NULL)
		ERROR("Environment variable VICTOR_ROOT was not found.\n Use the command:\n export VICTOR_ROOT=......", exception);
    string inputFile = getenv("VICTOR_ROOT");
  if (inputFile.length() < 3)
    ERROR("Environment variable VICTOR_ROOT was not found.", exception);
  
  inputFile += PARAM_FILE;
  
  ifstream input(inputFile.c_str());
  
  if (!input)
    ERROR("Could not read data file.", exception);
  
  vector <double>tmpA;
  vector<vector<double> >tmpB;
  for (int i = 0; i < 10; i++)
    tmpA.push_back(0);
  for (int i = 0; i < 10; i++)
    tmpB.push_back(tmpA);
  for (int i = 0; i < 21; i++)
    MaxPropensitiesDiscAngle.push_back(tmpB);
  
  int numCount;
  double propens;   
  string name;       
  int prephi, prepsi;

  input >> numCount;
  
  
  for (int i = 0; i < numCount; i++) 
    {
      input >> name >> prephi >> prepsi >> propens;
      AminoAcidCode r = aminoAcidThreeLetterTranslator(name);
      MaxPropensitiesDiscAngle[r][prephi][prepsi] = propens;
    }
  
  input.close();
  
  // calculate some constant values:
  
}

double MaxPropensitiesYPreAngle::GetPropensities (AminoAcidCode code, unsigned int prePhi , unsigned int prePsi )
{
  return MaxPropensitiesDiscAngle[code][prePhi][prePsi];
}



