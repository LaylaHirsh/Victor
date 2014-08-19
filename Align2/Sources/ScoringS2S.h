 

#ifndef __ScoringS2S_H__
#define __ScoringS2S_H__

#include <ScoringScheme.h>

namespace Biopool
{
/** @brief  Calculate scores for sequence to sequence alignment.
 * 
* @Description  
* @This 
 **/
class ScoringS2S : public ScoringScheme
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringS2S(SubMatrix *sub, AlignmentData *ad, Structure *str, double cSeq);

	/// Copy constructor.
	ScoringS2S(const ScoringS2S &orig);

	/// Destructor.
	virtual ~ScoringS2S();


// OPERATORS:

	/// Assignment operator.
	ScoringS2S& operator = (const ScoringS2S &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoring(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringS2S &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringS2S* newCopy();

	/// Reverse template sequence.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string seq1;    ///< Target sequence.
	string seq2;    ///< Template sequence.
	double cSeq;    ///< Coefficient for sequence alignment.

};

} // namespace

#endif
