 

#ifndef __Structure_H__
#define __Structure_H__

#include <SubMatrix.h>
#include <math.h>
#include <string>

namespace Biopool
{
/** @brief   Base class for structural scores.
 * 
* @Description  
* @This 
 **/
class Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Structure(SubMatrix *subStr) : subStr(subStr)
	{ }

	/// Copy constructor.
	Structure(const Structure &orig)
	{
		copy(orig);
	}

	/// Destructor.
	virtual ~Structure()
	{ }


// OPERATORS:

	/// Assignment operator.
	Structure& operator = (const Structure &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringStr(int i, int j) = 0;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Structure &orig);

	/// Construct a new "deep copy" of this object.
	virtual Structure* newCopy() = 0;

	/// Reverse template structural components.
	virtual void reverse()
	{ }


// ATTRIBUTES:

	SubMatrix *subStr;    ///< Structural substitution matrix.


protected:


private:

};

// -----------------------------------------------------------------------------
//                                  Structure
// -----------------------------------------------------------------------------

// OPERATORS:

inline Structure&
Structure::operator = (const Structure &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

inline void
Structure::copy(const Structure &orig)
{
	subStr = orig.subStr->newCopy();
}

} // namespace

#endif
