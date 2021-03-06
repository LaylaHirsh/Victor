/**
* @Class              SolvationPotential 
* @Project       Victor
*/

#ifndef _SOLVATIONPOTENTIAL_H_
#define _SOLVATIONPOTENTIAL_H_

// Includes:
#include <vector>
#include <Spacer.h>
#include <Potential.h>

// Global constants, typedefs, etc. (to avoid):
const double SOLVATION_CUTOFF_DISTANCE = 10.0;

namespace Biopool {
/** @brief class implements methods to get and set propensity
 * 
* @Description Includes methods that allow to  calculate solvation, energy and propensity.
* @This 
 * */
class SolvationPotential : public Potential{
public: 

  // CONSTRUCTORS/DESTRUCTOR:
  SolvationPotential(unsigned int _resol = 1);
  virtual ~SolvationPotential() { PRINT_NAME; } 

  // PREDICATES:
  virtual long double calculateEnergy(Spacer& sp) 
	{ return calculateSolvation(sp); }
  virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
				      unsigned int index2);
  virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp)
	{ return calculateSolvation(aa, sp); }
  virtual long double calculateEnergy(
    AminoAcid& resid, AminoAcidCode type, Spacer& sp) 
        { return calculateSolvation(resid, type, sp); }
  long double calculateSolvation(Spacer& sp);
  long double calculateSolvation(AminoAcid& aa, Spacer& sp, 
		  unsigned int start = 0, unsigned int end = 9999);
  long double calculateSolvation(AminoAcid& resid, AminoAcidCode type, 
                                 Spacer& sp);
  long double pReturnMaxPropensity(const AminoAcidCode type) const;
  long double pReturnMinPropensity(const AminoAcidCode type) const;

  // MODIFIERS:

  // OPERATORS:

  // HELPERS:

  static inline long double propCoeff() { return 0.582;}

protected:
  
private:

  // PREDICATES
  long double pGetPropensity(const AminoAcidCode type, 
                             unsigned int count) const;
  long double pGetMaxPropensity(const AminoAcidCode type) const;
  long double pGetMinPropensity(const AminoAcidCode type) const;

  // MODIFIERS
  void pConstructMaxPropensities();
  void pConstructMinPropensities();

  // ATTRIBUTES:
  vector<vector<int> > sum;
  
  unsigned int binResolution;

  /* maximum and minimum propensities extracted from the solvation database
     file. The i-th element of each vector refers to the i-th amino acid type
     declared in the AminoAcidCode enumeration (i=0,...,19). */
  //ALA=0,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR
  //    0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19  
  vector<long double> amino_max_propensities; 
  vector<long double> amino_min_propensities; 
  static unsigned int MAX_BINS;
  static string SOLV_PARAM_FILE;  


};

// ---------------------------------------------------------------------------
//                            SolvationPotential
// -----------------x-------------------x-------------------x-----------------
/**
 * @Description returns the maximum propensity for an amino acid type    
 * @param   the amino acid type(AminoAcidCode)
 * @return    the value of maximum propensity(long double)
 */
inline long double SolvationPotential::pReturnMaxPropensity(const AminoAcidCode type) const{
  return amino_max_propensities[type];
}
/**
 * @Description returns the minimum propensity for an amino acid type    
 * @param   the amino acid type(AminoAcidCode)
 * @return    the value of minimum propensity(long double)
 */
inline long double SolvationPotential::pReturnMinPropensity(const AminoAcidCode type) const{
  return amino_min_propensities[type];
}


} // namespace
#endif //_SOLVATIONPOTENTIAL_H_
