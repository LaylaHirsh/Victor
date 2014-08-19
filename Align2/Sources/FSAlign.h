 

#ifndef __FSAlign_H__
#define __FSAlign_H__

#include <Align.h>

namespace Biopool
{
/** @brief Implement free-shift "glocal" alignment.
 * 
* @Description  
* @This 
 **/
class FSAlign : public Align
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	FSAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss);

	/// Constructor with weighted alignment positions.
	FSAlign(AlignmentData *ad, GapFunction *gf, ScoringScheme *ss,
		const vector<unsigned int> &v1, const vector<unsigned int> &v2);

	/// Copy constructor.
	FSAlign(const FSAlign &orig);

	/// Destructor.
	virtual ~FSAlign();


// OPERATORS:

	/// Assignment operator.
	FSAlign& operator = (const FSAlign &orig);


// PREDICATES:

	/// Return two-element array containing an alignment with maximal score.
	virtual void getMultiMatch();


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const FSAlign &orig);

	/// Construct a new "deep copy" of this object.
	virtual FSAlign* newCopy();


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
