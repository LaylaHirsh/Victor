 

#ifndef __DotPFreq_H__
#define __DotPFreq_H__

#include <Profile.h>
#include <ScoringFunction.h>

namespace Biopool
{
/** @brief Calculate scores for profile to profile alignment using
*                  dot product method.
 * 
* @Description  Some explanations can be found in:
*
* 	                Mittelman D., Sadreyev R., Grishin N.
*                  Probabilistic scoring measures for profile-profile
*                  comparison yield more accurate short seed alignments.
*                  Bioinformatics. 2003 Aug 12;19(12):1531-9.
*                  PMID: 12912834 [PubMed - indexed for MEDLINE]
*
*                  Marti-Renom MA., Madhusudhan MS., Sali A.
*                  Alignment of protein sequences by their profiles.
*                  Protein Sci. 2004 Apr;13(4):1071-87.
*                  PMID: 15044736 [PubMed - indexed for MEDLINE]
* @This 
 **/
class DotPFreq : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	DotPFreq(Profile *pro1, Profile *pro2);

	/// Copy constructor.
	DotPFreq(const DotPFreq &orig);

	/// Destructor.
	virtual ~DotPFreq();


// OPERATORS:

	/// Assignment operator.
	DotPFreq& operator = (const DotPFreq &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const DotPFreq &orig);

	/// Construct a new "deep copy" of this object.
	virtual DotPFreq* newCopy();


protected:


private:

// ATTRIBUTES:

	Profile *pro1;    ///< Target profile.
	Profile *pro2;    ///< Template profile.

};

} // namespace

#endif
