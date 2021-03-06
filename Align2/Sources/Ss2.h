 

#ifndef __Ss2_H__
#define __Ss2_H__

#include <AlignmentData.h>
#include <Structure.h>
#include <Ss2Input.h>

namespace Biopool
{
/** @brief Calculate structural scores with info derived from PSI-PRED.  
 * 
* @Description  
* @This 
 **/
class Ss2 : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Ss2(SubMatrix *subStr, AlignmentData *ad, Ss2Input *psipred1,
		Ss2Input *psipred2, double cSs2);

	/// Copy constructor.
	Ss2(const Ss2 &orig);

	/// Destructor.
	virtual ~Ss2();


// OPERATORS:

	/// Assignment operator.
	Ss2& operator = (const Ss2 &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Ss2 &orig);

	/// Construct a new "deep copy" of this object.
	virtual Ss2* newCopy();

	/// Reverse template secondary structure.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string sec1;           ///< Target secondary structure.
	string sec2;           ///< Template secondary structure.
	Ss2Input *psipred1;    ///< Target PSI-PRED input file.
	Ss2Input *psipred2;    ///< Template PSI-PRED input file.
	double cSs2;           ///< Coefficient for PSI-PRED prediction.

};

} // namespace

#endif
