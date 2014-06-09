/**
* @Class              IRapdfPotential
* @Project         Victor
*/

#ifndef _IRAPDFPOTENTIAL_H_
#define _IRAPDFPOTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>
#include <RapdfPotential.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
/** @brief class implements the *interpolated* all-atom 
 * 
* @Description Includes methods that allow to manipulate the *interpolated* all-atom residue specific 
*probability discriminatory function from Samudrala & Moult (JMB 1998) .
* @This 
 * */
class IRapdfPotential : public RapdfPotential
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  IRapdfPotential();
  virtual ~IRapdfPotential() { PRINT_NAME; } 

  // PREDICATES:
  long double calculateEnergy(Spacer& sp);
  long double calculateEnergy(AminoAcid& aa, Spacer& sp);
  long double calculateEnergy(AminoAcid& aa, AminoAcid& aa2);
  long double calculateEnergy(Atom& at1, Atom& at2, string aaType, 
			      string aaType2);

  // MODIFIERS:

  // OPERATORS:

protected:

  // HELPERS:
  unsigned int pGetDistanceBinTwo(double distance);

private:

  // ATTRIBUTES:
  static string RAPDF_PARAM_FILE;

};

// ---------------------------------------------------------------------------
//                               IRapdfPotential
// -----------------x-------------------x-------------------x-----------------

inline long double 
IRapdfPotential::calculateEnergy(Atom& at1, Atom& at2, string aaType, 
string aaType2)
{
  double d = at1.distance(at2);
  double frac = d - static_cast<int>(d);

  if ((d >= 20.0) || (at1.getType() == "OXT") || (at2.getType() == "OXT")
      || ((at1.getType() == "CB") && (aaType == "GLY")) 
      || ((at2.getType() == "CB") && (aaType2 == "GLY")) )
    return 0.0;

  string tmp1 = threeLetter2OneLetter(aaType) + at1.getType();
  string tmp2 = threeLetter2OneLetter(aaType2) + at2.getType();

  return (1-frac) * prob[pGetDistanceBinOne(d)][pGetGroupBin(tmp1.c_str())][pGetGroupBin(tmp2.c_str())]
    +  frac * prob[pGetDistanceBinTwo(d)][pGetGroupBin(tmp1.c_str())][pGetGroupBin(tmp2.c_str())];
}

} // namespace
#endif //_IRAPDFPOTENTIAL_H_
