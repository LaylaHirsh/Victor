 

#ifndef __EDistance_H__
#define __EDistance_H__

#include <Profile.h>
#include <ScoringFunction.h>

namespace Biopool
{
/** @brief Calculate scores for profile to profile alignment using
*                  euclidean distance.
 * 
* @Description  
* @This 
 **/
class EDistance : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	EDistance(Profile *pro1, Profile *pro2);

	/// Constructor assigning offset.
	EDistance(Profile *pro1, Profile *pro2, double offset);

	/// Copy constructor.
	EDistance(const EDistance &orig);

	/// Destructor.
	virtual ~EDistance();


// OPERATORS:

	/// Assignment operator.
	EDistance& operator = (const EDistance &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);

	/// Return offset.
	virtual double getOffset();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const EDistance &orig);

	/// Construct a new "deep copy" of this object.
	virtual EDistance* newCopy();

	/// Set offset.
	virtual void setOffset(double off);


protected:


private:

// ATTRIBUTES:

	Profile *pro1;    ///< Target profile.
	Profile *pro2;    ///< Template profile.
	double offset;    ///< Offset.

};

// -----------------------------------------------------------------------------
//                                  EDistance
// -----------------------------------------------------------------------------

// PREDICATES:

inline double
EDistance::getOffset()
{
	return offset;
}


// MODIFIERS:

inline void
EDistance::setOffset(double off)
{
	offset = off;
}

} // namespace

#endif
