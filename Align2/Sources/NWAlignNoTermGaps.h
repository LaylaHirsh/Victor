 
#ifndef __NWAlignNoTermGaps_H__
#define __NWAlignNoTermGaps_H__

#include <Align.h>

namespace Biopool
{
/** @brief   Implement Needleman-Wunsch global alignment with no penalty
*                  for the terminal hangouts of the sequence.
 * 
* @Description  
* @This 
 **/
class NWAlignNoTermGaps : public Align
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	NWAlignNoTermGaps(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

	/// Constructor with weighted alignment positions.
	NWAlignNoTermGaps(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
		const vector<unsigned int> &v1, const vector<unsigned int> &v2);

	/// Copy constructor.
	NWAlignNoTermGaps(const NWAlignNoTermGaps &orig);

	/// Destructor.
	virtual ~NWAlignNoTermGaps();


// OPERATORS:

	/// Assignment operator.
	NWAlignNoTermGaps& operator = (const NWAlignNoTermGaps &orig);


// PREDICATES:

	/// Return two-element array containing an alignment with maximal score.
	virtual void getMultiMatch();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const NWAlignNoTermGaps &orig);

	/// Construct a new "deep copy" of this object.
	virtual NWAlignNoTermGaps* newCopy();


// HELPERS:

	/// Update/create matrix values.
	virtual void pCalculateMatrix(bool update = true);

	/// Update/create weighted matrix values.
	virtual void pCalculateMatrix(const vector<unsigned int> &v1,
		const vector<unsigned int> &v2, bool update = true);


protected:


private:

};

} // namespace

#endif
