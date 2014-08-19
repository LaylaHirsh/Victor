 

#ifndef __Pearson_H__
#define __Pearson_H__

#include <Profile.h>
#include <ScoringFunction.h>

namespace Biopool
{
/** @brief  Calculate scores for profile to profile alignment using
*                  Pearson's correlation coefficient. 
 * 
* @Description  Some explanations can be
*                  found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class Pearson : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Pearson(Profile *pro1, Profile *pro2);

	/// Copy constructor.
	Pearson(const Pearson &orig);

	/// Destructor.
	virtual ~Pearson();


// OPERATORS:

	/// Assignment operator.
	Pearson& operator = (const Pearson &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Pearson &orig);

	/// Construct a new "deep copy" of this object.
	virtual Pearson* newCopy();


protected:


private:

// ATTRIBUTES:

	Profile *pro1;    ///< Target profile.
	Profile *pro2;    ///< Template profile.
	double p1[20];    ///< Target background frequencies.
	double p2[20];    ///< Template background frequencies.

};

} // namespace

#endif
