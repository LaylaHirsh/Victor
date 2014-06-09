/**
 * @Class               EffectiveSolvationPotential
  * @Project        Victor
 * @Description 
*    This class implements a knowledge-based solvation with polar/hydrophobic 
*    information potential.
*/
#ifndef _EFFECTIVESOLVATIONPOTENTIAL_H_
#define _EFFECTIVESOLVATIONPOTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {

const double SOLVATION_CUTOFF_DISTANCE_EFFECTIVE = 10.0;
/** @brief class implements a knowledge-based solvation with polar/hydrophobic information potential
 * 
* @Description  
 * */
class EffectiveSolvationPotential : public Potential
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  EffectiveSolvationPotential();
  virtual ~EffectiveSolvationPotential() { PRINT_NAME; } 

  // PREDICATES:
  virtual long double calculateEnergy(Spacer& sp);
  virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
				      unsigned int index2);
  virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp);

  // MODIFIERS:

  // OPERATORS:

protected:

  // HELPERS:
  bool isPolar(AminoAcid& aa);
  double pCalcFracBuried(unsigned int index, Spacer& sp);

private:

  // ATTRIBUTES:

  vector<double> solvCoeff;

};

// ---------------------------------------------------------------------------
//                          EffectiveSolvationPotential
// -----------------x-------------------x-------------------x-----------------



/**
 * @Description Verifies if the amino acid is a polar one
 * @param reference of the aa(aminoAcid&)
 * @return  result of the validation(bool)
 */

inline bool EffectiveSolvationPotential::isPolar(AminoAcid& aa){
  AminoAcidCode c = static_cast<AminoAcidCode>( aa.getCode() );
  if ( (c == ARG) ||
       (c == ASN) ||
       (c == ASP) ||
       (c == GLN) ||
       (c == GLU) ||
       (c == LYS) ||
       (c == PRO) )
    return true;

  return false;
}


} // namespace
#endif //_EFFECTIVESOLVATIONPOTENTIAL_H_
