 
#ifndef __NWAlign_H__
#define __NWAlign_H__

#include <Align.h>

namespace Biopool
{
/** @brief  Implement Needleman-Wunsch global alignment.
 * 
* @Description  
* @This 
 **/
class NWAlign : public Align
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	NWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

	/// Constructor with weighted alignment positions.
	NWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
		const vector<unsigned int> &v1, const vector<unsigned int> &v2);

	/// Copy constructor.
	NWAlign(const NWAlign &orig);

	/// Destructor.
	virtual ~NWAlign();


// OPERATORS:

	/// Assignment operator.
	NWAlign& operator = (const NWAlign &orig);


// PREDICATES:

	/// Return two-element array containing an alignment with maximal score.
	virtual void getMultiMatch();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const NWAlign &orig);

	/// Construct a new "deep copy" of this object.
	virtual NWAlign* newCopy();


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
