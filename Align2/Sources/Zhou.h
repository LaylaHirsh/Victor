 

#ifndef __Zhou_H__
#define __Zhou_H__

#include <Profile.h>
#include <PssmInput.h>
#include <ScoringFunction.h>

namespace Biopool
{

/** @brief   Calculate scores for profile to profile alignment using
//                  Zhou-Zhou method.
 * 
* @Description   Some explanations can be found in:
*
*                  Zhou H., Zhou Y.
*                  Single-body residue-level knowledge-based energy score
*                  combined with sequence-profile and secondary structure
*                  information for fold recognition.
*                  Proteins. 2004 Jun 1;55(4):1005-13.
*                  PMID: 15146497 [PubMed - indexed for MEDLINE]
* @This 
 **/
class Zhou : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Zhou(Profile *pro1, PssmInput *pssm2);

	/// Copy constructor.
	Zhou(const Zhou &orig);

	/// Destructor.
	virtual ~Zhou();


// OPERATORS:

	/// Assignment operator.
	Zhou& operator = (const Zhou &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Zhou &orig);

	/// Construct a new "deep copy" of this object.
	virtual Zhou* newCopy();


protected:


private:

// ATTRIBUTES:

	Profile *pro1;       ///< Target profile.
	PssmInput *pssm2;    ///< Template PSSM.

};

} // namespace

#endif
