

              /*************************************************
               * sub is BLOSUM62 matrix standard log-odds form *
               *************************************************/

#ifndef __CrossProduct_H__
#define __CrossProduct_H__

#include <Profile.h>
#include <ScoringFunction.h>
#include <SubMatrix.h>

namespace Biopool
{
/** @brief  Calculate scores for profile to profile alignment using
*                  sum of pairs method. 
 * 
* @Description  Some explanations can be found in:
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class CrossProduct : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	CrossProduct(SubMatrix *sub, Profile *pro1, Profile *pro2);

	/// Copy constructor.
	CrossProduct(const CrossProduct &orig);

	/// Destructor.
	virtual ~CrossProduct();


// OPERATORS:

	/// Assignment operator.
	CrossProduct& operator = (const CrossProduct &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const CrossProduct &orig);

	/// Construct a new "deep copy" of this object.
	virtual CrossProduct* newCopy();


protected:


private:

// ATTRIBUTES:

	SubMatrix *sub;    ///< Substitution matrix.
	Profile *pro1;     ///< Target profile.
	Profile *pro2;     ///< Template profile.

};

} // namespace

#endif
