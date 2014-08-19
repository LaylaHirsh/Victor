 

#ifndef __ScoringP2S_H__
#define __ScoringP2S_H__

#include <Profile.h>
#include <ScoringScheme.h>

namespace Biopool
{
/** @brief   Calculate scores for profile to sequence alignment.
 * 
* @Description  
* @This 
 **/
class ScoringP2S : public ScoringScheme
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringP2S(SubMatrix *sub, AlignmentData *ad, Structure *str, Profile *pro,
		double cSeq);

	/// Copy constructor.
	ScoringP2S(const ScoringP2S &orig);

	/// Destructor.
	virtual ~ScoringP2S();


// OPERATORS:

	/// Assignment operator.
	ScoringP2S& operator = (const ScoringP2S &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoring(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringP2S &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringP2S* newCopy();

	/// Reverse template sequence.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string seq1;     ///< Target sequence.
	string seq2; 	 ///< Template sequence.
	Profile *pro;    ///< Target profile.
	double cSeq;     ///< Coefficient for sequence alignment.

};

} // namespace

#endif
