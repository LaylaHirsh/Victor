 /**
 * @Class               EnergyFeatures 
 * @Project        Victor
 * @Description 
*    Interface/wrapper for energy feature calculation, e.g. in FRST2.
*/
#ifndef _ENERGYFEATURE_H_
#define _ENERGYFEATURE_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <RapdfPotential.h>
#include <SolvationPotential.h>
#include <PhiPsi.h>
#include <PhiPsiOmegaChi1Chi2.h>


// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief  Interface/wrapper for energy feature calculation, e.g. in FRST2.
 * 
* @Description 
* @This 
 * */
class EnergyFeatures{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  EnergyFeatures();
  virtual ~EnergyFeatures() { PRINT_NAME; } 

  // PREDICATES:
  static void showFeatures();
  static double calculateBackboneHydrogenBonds(Spacer& sp);
  static vector<double> calculateAaComposition(Spacer& sp);
  static vector<double> calculateSecondaryComposition(Spacer& sp);
  static vector<double> calculateMesoStateComposition(Spacer& sp);
  static vector<double> calculateChainBreaks(Spacer& sp);
  static double calculateClashes(Spacer& sp);

  /// Main function wrapping up feature calculation for FRST2
  vector<double> calculateFeatures(Spacer& sp);

  // MODIFIERS:

  // OPERATORS:

protected:

private:

  // HELPERS:

  // ATTRIBUTES:

  SolvationPotential solv;
  RapdfPotential rapdf;
  TorsionPotential* tors;
  TorsionPotential* tors2;

};

// ---------------------------------------------------------------------------
//                                EnergyFeatures
// -----------------x-------------------x-------------------x-----------------


} // namespace
#endif //_ENERGYFEATURE_H_
