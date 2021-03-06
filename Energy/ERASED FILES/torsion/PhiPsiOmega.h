/**
*@Class               PhiPsiOmega
*@Modified        Andrea Bazzoli
*@Project        Victor
*/

#ifndef _PHIPSIOMEGA_H_
#define _PHIPSIOMEGA_H_


// Includes:
#include <TorsionPotential.h>
#include <vector>
#include <Spacer.h>
#include <AminoAcidCode.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
  /** @brief class manages the angle qualities and the energy 
 * 
* @Description  This class implements a simple torsion potential based on the statistical preference of aminoacid types phi , psi and omega angles.
* @This 
 * */
	class PhiPsiOmega : public TorsionPotential
{
public: 
	
	// CONSTRUCTORS/DESTRUCTOR:
	PhiPsiOmega ( int SET_ARC1 = 10,
				  string knownledge = "data/tor.par" );//default knownledge TOP500
	virtual ~PhiPsiOmega() { pResetData(); } 
	
	// PREDICATES:
	virtual long double calculateEnergy(Spacer& sp);
	virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
										unsigned int index2);
	virtual long double calculateEnergy(AminoAcid& aa);
	virtual long double calculateEnergy(AminoAcid& diheds, AminoAcidCode code);
	virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp)
    { return calculateEnergy(aa); }
	virtual double pReturnMaxPropensities(int amino);
	virtual double pReturnMinPropensities(int amino);
	virtual double pReturnMaxPropensitiesPreAngle(int amino, int prephi, int prepsi)
	{ ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
	virtual int sGetPropBin2(double p)
	{ ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.",exception)}
   	virtual vector< vector<ANGLES> >* getOrderedEnergyTable();
	
	string getLabel() {return "phi-psi-omega";}
	
	// MODIFIERS:
	virtual int setRange_Omega(int n);
	virtual void setArcStep(int n);
	
	// OPERATORS:
	
protected:
		
		// HELPERS:
		virtual void pConstructData();
	virtual void pResetData();
	virtual double pGetMaxPropensities( int amino );
	virtual double pGetMinPropensities( int amino ); 
	//not implemented in this class. No pre-angle considered.
	virtual double pGetMaxPropensities( int amino, int prephi, int prepsi)
	{ ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
	void sAddProp(int code, int x, int y, int z);
	int sGetPropOmegaBin(double p);
	int sGetPropBin(double p);
	virtual void pConstructMaxPropensities();
	virtual void pConstructMinPropensities();
	virtual double getOmegaAngle( int prop, long RANGE_OMEGA );
	virtual double getPhiPsiAngle( int prop, long SIZE_TABLE, int ARC_STEP );
	
	
	//private:
public:
		
		// ATTRIBUTES:
		int ARC_STEP;  // important: must be a divisior of 360 !!!! 
	string TOR_PARAM_FILE;   // File with prop torsion angles
	int SIZE_OF_TABLE; // "granularity" props
	int RANGE_OMEGA;
	int amino_count[AminoAcid_CODE_SIZE];//number of amino. 
		vector<vector<vector<vector<int>* >* >* > propensities;// the propensities table.
			vector<vector<vector<int>* >* > all_propensities;// the sum of propropensities table.
				double total;//total numer of ammino considered.
					vector<double> amino_max_propensities; // vectors with max and min amino propensities
					vector<double> amino_min_propensities; // according to knowledge. 
};

// ---------------------------------------------------------------------------
//                            TorsionPotential
// -----------------x-------------------x-------------------x-----------------
} // namespace
#endif// _FULLTORSIONPOTENTIAL_H_












