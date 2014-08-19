 

#ifndef __ScoringScheme_H__
#define __ScoringScheme_H__

#include <AlignmentData.h>
#include <Structure.h>
#include <SubMatrix.h>
#include <math.h>
#include <string>

namespace Biopool
{
/** @brief  Base class for scoring schemes.
 * 
* @Description  
* @This 
 **/
class ScoringScheme
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringScheme(SubMatrix *sub, AlignmentData *ad, Structure *str);

	/// Copy constructor.
	ScoringScheme(const ScoringScheme &orig);

	/// Destructor.
	virtual ~ScoringScheme();


// OPERATORS:

	/// Assignment operator.
	ScoringScheme& operator = (const ScoringScheme &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoring(int i, int j) = 0;

	/// Check if s consists only of characters defined in sub.getResidues.
	virtual bool checkSequence(const string &s) const;


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringScheme &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringScheme* newCopy() = 0;

	/// Reverse template components (sequence and/or profile).
	virtual void reverse();


// ATTRIBUTES:

	SubMatrix *sub;       ///< Substitution matrix.
	AlignmentData *ad;    ///< Pointer to AlignmentData.
	Structure *str;       ///< Pointer to Structure.


protected:


private:

};

} // namespace

#endif
