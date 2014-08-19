 

#ifndef __GapFunction_H__
#define __GapFunction_H__

#include <Debug.h>
#include <iostream>

namespace Biopool
{
/** @brief Base class for gap functions.
 * 
* @Description  
* @This 
 **/
class GapFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	GapFunction()
	{ }

	/// Copy constructor.
	GapFunction(const GapFunction &orig)
	{
		copy(orig);
	}

	/// Destructor.
	virtual ~GapFunction()
	{ }


// OPERATORS:

	/// Assignment operator.
	GapFunction& operator = (const GapFunction &orig);


// PREDICATES:

	/// Return open gap penalty for template position p.
	virtual double getOpenPenalty(int p) = 0;

	/// Return extension gap penalty for template position p.
	virtual double getExtensionPenalty(int p) = 0;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const GapFunction &orig);

	/// Construct a new "deep copy" of this object.
	virtual GapFunction* newCopy() = 0;

	/// Set open gap penalty.
	virtual void setOpenPenalty(double pen) = 0;

	/// Set extension gap penalty.
	virtual void setExtensionPenalty(double pen) = 0;


protected:


private:

};

// -----------------------------------------------------------------------------
//                                 GapFunction
// -----------------------------------------------------------------------------

// OPERATORS:

inline GapFunction&
GapFunction::operator = (const GapFunction &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// MODIFIERS:

inline void
GapFunction::copy(const GapFunction &orig)
{ }

} // namespace

#endif
