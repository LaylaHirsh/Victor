 
#ifndef __Sec_H__
#define __Sec_H__

#include <AlignmentData.h>
#include <Structure.h>

namespace Biopool
{
/** @brief   Calculate structural scores with info derived from secondary
*                  structure.
 * 
* @Description  
* @This 
 **/
class Sec : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Sec(SubMatrix *subStr, AlignmentData *ad, double cSec);

	/// Copy constructor.
	Sec(const Sec &orig);

	/// Destructor.
	virtual ~Sec();


// OPERATORS:

	/// Assignment operator.
	Sec& operator = (const Sec &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Sec &orig);

	/// Construct a new "deep copy" of this object.
	virtual Sec* newCopy();

	/// Reverse template secondary structure.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string sec1;           ///< Target secondary structure.
	string sec2;           ///< Template secondary structure.
	double cSec;           ///< Coefficient for secondary structure alignment.

};

} // namespace

#endif
