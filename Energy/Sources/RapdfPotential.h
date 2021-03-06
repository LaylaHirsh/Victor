/**
* @Class               RapdfPotential 
* @Project       Victor 
**/

#ifndef _RAPDFPOTENTIAL_H_
#define _RAPDFPOTENTIAL_H_
// Includes:
#include <vector>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):
const unsigned int MAX_BINS = 18;
const unsigned int MAX_TYPES = 168;

namespace Biopool {
/** @brief class implements the all-atom residue.
 * 
* @Description Includes methods that allow to manipulate the all-atom residue specific probability 
//    discriminatory function from Samudrala & Moult (JMB 1998).
 * */
class RapdfPotential : public Potential
{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  RapdfPotential();
  virtual ~RapdfPotential() { PRINT_NAME; } 

  // PREDICATES:
  virtual long double calculateEnergy(Spacer& sp);
  virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
				      unsigned int index2);
  virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp);
  virtual long double calculateEnergy(AminoAcid& aa, AminoAcid& aa2);
  virtual long double calculateEnergy(Atom& at1, Atom& at2, string aaType, 
			      string aaType2);
 
  // MODIFIERS: 
  // OPERATORS:

protected:

  // HELPERS:
  unsigned int pGetDistanceBinOne(double distance);
  unsigned int pGetGroupBin(const char* group_name);
  
  // ATTRIBUTES:
  double prob[MAX_BINS][MAX_TYPES][MAX_TYPES];

public:
    //static string RAPDF_PARAM_FILE; //error at runtime in a 64 bit architecture, uncomment this if you are in a no 64bits SO
    string path;				       
private:

};

// ---------------------------------------------------------------------------
//                               RapdfPotential
// -----------------x-------------------x-------------------x-----------------


/**
 * @Description calculates the energy between two atoms 
 * @param   the references to the atoms (Atom&,Atom&)the amino acid types(string,string)
 * @return    the value of maximum propensity(long double)
 */
    inline long double RapdfPotential::calculateEnergy(Atom& at1, Atom& at2, string aaType, 
    string aaType2){
      double d = at1.distance(at2);
      if ((d >= 20.0) || (at1.getType() == "OXT") || (at2.getType() == "OXT")
          || ((at1.getType() == "CB") && (aaType == "GLY")) 
          || ((at2.getType() == "CB") && (aaType2 == "GLY")) )
        return 0.0;
      string tmp1 = threeLetter2OneLetter(aaType) + at1.getType();
      string tmp2 = threeLetter2OneLetter(aaType2) + at2.getType();
      unsigned int dist = pGetDistanceBinOne(d);
      unsigned int grp1 = pGetGroupBin(tmp1.c_str());
      unsigned int grp2 = pGetGroupBin(tmp2.c_str());
      if (dist + grp1 + grp2 < 999)   
        return prob[dist][grp1][grp2];
      else // ignore errors
        return 0;
    }

} // namespace
#endif //_RAPDFPOTENTIAL_H_
