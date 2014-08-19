 

#ifndef __ScoringP2P_H__
#define __ScoringP2P_H__

#include <Profile.h>
#include <ScoringFunction.h>
#include <ScoringScheme.h>

namespace Biopool
{
/** @brief  Calculate scores for profile to profile alignment.
 * 
* @Description  
* @This 
 **/
class ScoringP2P : public ScoringScheme
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ScoringP2P(SubMatrix *sub, AlignmentData *ad, Structure *str, Profile *pro1,
		Profile *pro2, ScoringFunction *fun, double cSeq);

	/// Copy constructor.
	ScoringP2P(const ScoringP2P &orig);

	/// Destructor.
	virtual ~ScoringP2P();


// OPERATORS:

	/// Assignment operator.
	ScoringP2P& operator = (const ScoringP2P &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoring(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ScoringP2P &orig);

	/// Construct a new "deep copy" of this object.
	virtual ScoringP2P* newCopy();

	/// Reverse template sequence and profile.
	virtual void reverse();


protected:


private:

// ATTRIBUTES:

	string seq1;             ///< Target sequence.
	string seq2;             ///< Template sequence.
	Profile *pro1;           ///< Target profile.
	Profile *pro2;           ///< Template profile.
	ScoringFunction *fun;    ///< Scoring function.
	double cSeq;             ///< Coefficient for sequence alignment.

};

} // namespace

#endif
