 

#ifndef __Prof_H__
#define __Prof_H__

#include <ProfInput.h>
#include <Structure.h>

namespace Biopool
{
/** @brief   Calculate structural scores with info derived from PHD.
 * 
* @Description  
* @This 
 **/
class Prof : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Prof(SubMatrix *subStr, ProfInput *phd1, ProfInput *phd2, double cPrf);

	/// Copy constructor.
	Prof(const Prof &orig);

	/// Destructor.
	virtual ~Prof();


// OPERATORS:

	/// Assignment operator.
	Prof& operator = (const Prof &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Prof &orig);

	/// Construct a new "deep copy" of this object.
	virtual Prof* newCopy();


protected:


private:

// ATTRIBUTES:

	ProfInput *phd1;    ///< Target PHD input file.
	ProfInput *phd2;    ///< Template PHD input file.
	double cPrf;        ///< Coefficient for PHD prediction.

};

} // namespace

#endif
