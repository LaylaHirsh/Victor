  /**
 * @Class             Potential 
 * @Description 
*    Base class for energy calculation
*/
#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_


// Includes:
#include <vector>
#include <Spacer.h>

// Global constants, typedefs, etc. (to avoid):

namespace Biopool {
    /** @brief class contains the angle structure and a vector with the ordered energies 
 * 
* @Description Includes methods that allow to obtain propensity. 
 * */
	class Potential
{
public: 
	struct ANGLES
	{
		double phi;
		double psi;
		double omega;
		double energy;
	};
	
	// CONSTRUCTORS/DESTRUCTOR:
	Potential() { }
	Potential(string inputFile) { }
	virtual ~Potential() { PRINT_NAME; } 
	
	// PREDICATES:
	virtual long double calculateEnergy(Spacer& sp) = 0;
	virtual long double calculateEnergy(Spacer& sp, unsigned int index1, 
										unsigned int index2) = 0;
	virtual long double calculateEnergy(AminoAcid& aa, Spacer& sp) = 0;
	virtual long double calculateEnergy(AminoAcid& resid, AminoAcidCode type, Spacer& sp)
		{ ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception)}
	
	virtual long double pReturnMaxPropensity(const AminoAcidCode type) const 
		{ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);}
	virtual long double pReturnMinPropensity(const AminoAcidCode type) const 
		{ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception);}
	
	virtual vector< vector<ANGLES> >* orderedEnergy()
		{ERROR("ERROR. NOT IMPLEMENTED FOR THIS CLASS.", exception); }
	// MODIFIERS:
	
	// OPERATORS:
	
protected:
		
		// HELPERS:
		
		// ATTRIBUTES:
		
private:
		
};

// ---------------------------------------------------------------------------
//                                  Potential
// -----------------x-------------------x-------------------x-----------------
    /** @brief Energy operator
 * */
struct EnergyGreater : public binary_function<const Potential::ANGLES &,const Potential::ANGLES &,bool>
{
	bool operator()( const Potential::ANGLES &ref1, const Potential::ANGLES &ref2 )
	{
		return ref1.energy>ref2.energy;
	}
};

} // namespace
#endif //_POTENTIAL_H_
