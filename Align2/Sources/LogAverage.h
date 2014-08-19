 

              /*************************************************
               * sub is BLOSUM62 matrix standard log-odds form *
               *************************************************/

#ifndef __LogAverage_H__
#define __LogAverage_H__

#include <Profile.h>
#include <ScoringFunction.h>
#include <SubMatrix.h>

namespace Biopool
{
/** @brief Calculate scores for profile to profile alignment using
*                  sum of pairs method. 
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class LogAverage : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	LogAverage(SubMatrix *sub, Profile *pro1, Profile *pro2);

	/// Copy constructor.
	LogAverage(const LogAverage &orig);

	/// Destructor.
	virtual ~LogAverage();


// OPERATORS:

	/// Assignment operator.
	LogAverage& operator = (const LogAverage &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const LogAverage &orig);

	/// Construct a new "deep copy" of this object.
	virtual LogAverage* newCopy();


protected:


private:

// ATTRIBUTES:

	SubMatrix *sub;    ///< Substitution matrix.
	Profile *pro1;     ///< Target profile.
	Profile *pro2;     ///< Template profile.

};

} // namespace

#endif
