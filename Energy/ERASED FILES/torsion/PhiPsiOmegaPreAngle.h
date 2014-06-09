/**
*@Class               PhiPsiOmegaPreAngle
*@Base      TorsionPotential
*@Project         Victor
*/

#ifndef _PHIPSIOMEGAPREANGLE_H_
#define _PHIPSIOMEGAPREANGLE_H_


// Includes:
#include <TorsionPotential.h>
#include <vector>
#include <Spacer.h>
#include <AminoAcidCode.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
  /** @brief class manages the angle qualities and the energy 
 * 
* @Description This class implements a simple torsion potential based on the statistical preference of aminoacid types for  phi, psi, omega, prephi and pre psi angles.
* @This 
 * */
  class PhiPsiOmegaPreAngle : public TorsionPotential
  {
  public: 
    
    // CONSTRUCTORS/DESTRUCTOR:
    PhiPsiOmegaPreAngle ( int SET_ARC1 = 10, int SET_ARC2 = 40, 
			  string knownledge = "data/tor.par");
    virtual ~PhiPsiOmegaPreAngle( ) { pResetData(); } 
    
    // PREDICATES:
    virtual long double calculateEnergy(Spacer& sp);
    virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
					unsigned int index2);
    virtual long double calculateEnergy(AminoAcid& aa);
    virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp)
    { return calculateEnergy( aa ); }
    virtual long double calculateEnergy(AminoAcid& diheds, AminoAcidCode code) 
    {ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
  
    virtual double pReturnMaxPropensities(int amino);
    virtual double pReturnMaxPropensitiesPreAngle(int amino, int prephi, int prepsi);
    virtual int sGetPropBin2(double p);
    // MODIFIERS:
    virtual void setArcStep(int n);
    virtual int setRange_Omega(int n);

    // OPERATORS:
    
  protected:
    
    // HELPERS:
    virtual void pConstructData();
    virtual void pResetData();
    virtual double pGetMaxPropensities(int amino);
    virtual double pGetMaxPropensities(int amino, int prephi, int prepsi);
    virtual void sAddProp(int code, int x, int y, int z, int m, int n);
    virtual int sGetPropBin (double p);
    virtual int sGetPropOmegaBin(double p);
    virtual void pConstructMaxPropensities();
  private:
    
      // ATTRIBUTES:
    string TOR_PARAM_FILE;   // File with prop torsion angles
    int RANGE_OMEGA;
    int ARC_STEP;  // important: must be a divisior of 360 !!!!
    int ARC_STEP2;// the Arcstep for the pre-angles
    int SIZE_OF_TABLE; // "granularity" props
    int SIZE_OF_TABLE2;//"granularity" props for pre-angles
    int amino_count[AminoAcid_CODE_SIZE]; // total number of entries for all amino acids
    vector<vector<vector<vector<vector<vector<int>* >* >* >* >* > propensities;
    vector<vector<vector<vector<vector<int>* >* >* >* > all_propensities;
    double total;
    vector<double> amino_max_propensities;//vector with max amino propensities
                                        // according to knowledge.
    vector<vector<vector<double>* >* > amino_max_propensities_pre_angle;
    
  };
  
  // ---------------------------------------------------------------------------
  //                            PhiPsiOmegaPreAngle
  // -----------------x-------------------x-------------------x-----------------
} // namespace
#endif// _PHIPSIOMEGAPREANGLE_H_


