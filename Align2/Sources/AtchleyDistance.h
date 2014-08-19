

#ifndef __AtchleyDistance_H__
#define __AtchleyDistance_H__

#include <Profile.h>
#include <ScoringFunction.h>
#include <iostream>

namespace Biopool
{
/** @brief  Calculate scores for profile to profile alignment using
*                  sequence metric factor.
 * 
* @Description   Some explanations can be found in:
*
*                  William R. Atchley, Jieping Zhao, Andrew D. Fernandes, Tanja Druke
*                  Solving the protein sequence metric problem.
*                  Edited by Walter M. Fitch, University of California, Irvine, CA,
*                  and approved March 22, 2005 (received for review December 14, 2004).
* @This 
 **/
class AtchleyDistance : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	AtchleyDistance(Profile *pro1, Profile *pro2);

	/// Constructor assigning offset.
	AtchleyDistance(Profile *pro1, Profile *pro2, double offset);

	/// Copy constructor.
	AtchleyDistance(const AtchleyDistance &orig);

	/// Destructor.
	virtual ~AtchleyDistance();


// OPERATORS:

	/// Assignment operator.
	AtchleyDistance& operator = (const AtchleyDistance &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);

	/// Return offset.
	virtual double getOffset();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const AtchleyDistance &orig);

	/// Construct a new "deep copy" of this object.
	virtual AtchleyDistance* newCopy();

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
//                               AtchleyDistance
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
AtchleyDistance::getOffset()
{
	return offset;
}


// MODIFIERS:

inline void
AtchleyDistance::setOffset(double off)
{
	offset = off;
}

} // namespace

#endif
