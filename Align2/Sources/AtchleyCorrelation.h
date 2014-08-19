

#ifndef __AtchleyCorrelation_H__
#define __AtchleyCorrelation_H__

#include <Profile.h>
#include <ScoringFunction.h>
#include <iostream>

namespace Biopool
{
/** @brief  Calculate scores for profile to profile alignment using
*                  sequence metric factor. 
 * 
* @Description  Some explanations can be found in:
*
*                  William R. Atchley, Jieping Zhao, Andrew D. Fernandes, Tanja Druke
*                  Solving the protein sequence metric problem.
*                  Edited by Walter M. Fitch, University of California, Irvine, CA,
*                  and approved March 22, 2005 (received for review December 14, 2004).
* @This 
 **/
class AtchleyCorrelation : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	AtchleyCorrelation(Profile *pro1, Profile *pro2);

	/// Constructor assigning offset.
	AtchleyCorrelation(Profile *pro1, Profile *pro2, double offset);

	/// Copy constructor.
	AtchleyCorrelation(const AtchleyCorrelation &orig);

	/// Destructor.
	virtual ~AtchleyCorrelation();


// OPERATORS:

	/// Assignment operator.
	AtchleyCorrelation& operator = (const AtchleyCorrelation &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);

	/// Return offset.
	virtual double getOffset();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const AtchleyCorrelation &orig);

	/// Construct a new "deep copy" of this object.
	virtual AtchleyCorrelation* newCopy();

	/// Set offset.
	virtual void setOffset(double off);


// HELPERS:

	/// Helper function used to load Atchley metric factor.
	virtual void pLoadFactor();


protected:


private:

// ATTRIBUTES:

	Profile *pro1;           ///< Target profile.
	Profile *pro2;           ///< Template profile.
	double offset;           ///< Offset.
	double factor[20][5];    ///< Atchley metric factor.

};

// -----------------------------------------------------------------------------
//                              AtchleyCorrelation
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
AtchleyCorrelation::getOffset()
{
	return offset;
}


// MODIFIERS:

inline void
AtchleyCorrelation::setOffset(double off)
{
	offset = off;
}

} // namespace

#endif
