 

#ifndef __SWAlign_H__
#define __SWAlign_H__

#include <Align.h>

namespace Biopool
{
/** @brief  Implement Smith-Waterman local alignment.
 * 
* @Description  
* @This 
 **/
class SWAlign : public Align
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

	/// Constructor with weighted alignment positions.
	SWAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
		const vector<unsigned int> &v1, const vector<unsigned int> &v2);

	/// Copy constructor.
	SWAlign(const SWAlign &orig);

	/// Destructor.
	virtual ~SWAlign();


// OPERATORS:

	/// Assignment operator.
	SWAlign& operator = (const SWAlign &orig);


// PREDICATES:

	/// Return two-element array containing an alignment with maximal score.
	virtual void getMultiMatch();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SWAlign &orig);

	/// Construct a new "deep copy" of this object.
	virtual SWAlign* newCopy();


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
