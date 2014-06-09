 /**
 * @Class               PolarSolvationPotential 
 * @Project        Victor
*/
#ifndef _POLARSOLVATIONPOTENTIAL_H_
#define _POLARSOLVATIONPOTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

const double SOLVATION_CUTOFF_DISTANCE_POLAR = 7.0;
/** @brief This class implements a knowledge-based solvation with polar/hydrophobic   information potential.
 * 
 * */
class PolarSolvationPotential : public Potential
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  PolarSolvationPotential();
  virtual ~PolarSolvationPotential() { PRINT_NAME; } 

  // PREDICATES:
  virtual long double calculateEnergy(Spacer& sp) 
	{ return calculateSolvation(sp); }
  virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
				      unsigned int index2);
  virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp)
	{ return calculateSolvation(aa, sp); }
  long double calculateSolvation(Spacer& sp);
  long double calculateSolvation(AminoAcid& aa, Spacer& sp, 
				 unsigned int start = 0, 
				 unsigned int end = 9999);

  // MODIFIERS:

  // OPERATORS:

protected:

  // HELPERS:

private:

  // ATTRIBUTES:
  vector<vector<vector<int> > > sumPolar;

  static unsigned int MAX_BINS;
  static unsigned int BIN_POLAR;
  static string SOLV_PARAM_FILE;
};

} // namespace
#endif //_POLARSOLVATIONPOTENTIAL_H_
