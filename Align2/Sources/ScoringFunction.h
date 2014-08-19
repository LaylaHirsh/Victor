 

#ifndef __ScoringFunction_H__
#define __ScoringFunction_H__

#include <math.h>
#include <string>

namespace Biopool
{
/** @brief    Base class for scoring functions.
 * 
* @Description  
* @This 
 **/
class ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringFunction()
	{ }

	/// Copy constructor.
	ScoringFunction(const ScoringFunction &orig)
	{
		copy(orig);
	}

	/// Destructor.
	virtual ~ScoringFunction()
	{ }


// OPERATORS:

	/// Assignment operator.
	ScoringFunction& operator = (const ScoringFunction &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j) = 0;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringFunction &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringFunction* newCopy() = 0;


protected:


private:

};

// -----------------------------------------------------------------------------
//                               ScoringFunction
// -----------------------------------------------------------------------------

// OPERATORS:

inline ScoringFunction&
ScoringFunction::operator = (const ScoringFunction &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

inline void
ScoringFunction::copy(const ScoringFunction &orig)
{ }

} // namespace

#endif
