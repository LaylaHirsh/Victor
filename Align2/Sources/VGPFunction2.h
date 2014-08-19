 

#ifndef __VGPFunction2_H__
#define __VGPFunction2_H__

#include <Alignment.h>
#include <GapFunction.h>
#include <math.h>
#include <string>
#include <vector>

namespace Biopool
{
    /** @brief   Implement VGP (Variable Gap Penalty) function.
 * 
* @Description   Some
*                  explanations can be found in:
*
*                  Madhusudhan MS., Marti-Renom MA., Sanchez R., Sali A.
*                  Variable gap penalty for protein sequence-structure alignment.
*                  Department of Biopharmaceutical Sciences and Pharmaceutical
*                  Chemistry, University of California at San Francisco, 94143, USA.
*                  PMID: 16423846 [PubMed - indexed for MEDLINE]
* @This 
 **/

class VGPFunction2 : public GapFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	VGPFunction2(string secFileName);

	/// Constructor assigning o, e and weights.
	VGPFunction2(string secFileName, double o, double e, unsigned int extType,
		double wH, double wS);

	/// Copy constructor.
	VGPFunction2(const VGPFunction2 &orig);

	/// Destructor.
	virtual ~VGPFunction2();


// OPERATORS:

	/// Assignment operator.
	VGPFunction2& operator = (const VGPFunction2 &orig);


// PREDICATES:

	/// Return open gap penalty for template position p.
	virtual double getOpenPenalty(int p);

	/// Return extension gap penalty for template position p.
	virtual double getExtensionPenalty(int p);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const VGPFunction2 &orig);

	/// Construct a new "deep copy" of this object.
	virtual VGPFunction2* newCopy();

	/// Set open gap penalty.
	virtual void setOpenPenalty(double pen);

	/// Set extension gap penalty.
	virtual void setExtensionPenalty(double pen);


// HELPERS:

	/// Extract structural infos from template secondary structure.
	void pExtractSecInfo(string secFileName);


protected:


private:

// ATTRIBUTES:

	double o;                     ///< Open gap penalty.
	double e;                     ///< Extension gap penalty.
	unsigned int extType;         ///< Extension gap type.
	unsigned int extCounter;      ///< Extension gap counter.
	vector<double> hContent;      ///< Template helical content.
	vector<double> sContent;      ///< Template strand content.
	double wH;                    ///< Weight for helical content.
	double wS;                    ///< Weight for strand content.

};

// -----------------------------------------------------------------------------
//                                 VGPFunction
// -----------------------------------------------------------------------------

// MODIFIERS:

inline void
VGPFunction2::setOpenPenalty(double pen)
{
	o = pen;
}


inline void
VGPFunction2::setExtensionPenalty(double pen)
{
	e = pen;
}

} // namespace

#endif
