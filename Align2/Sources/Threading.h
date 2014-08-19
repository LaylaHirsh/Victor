 

#ifndef __Threading_H__
#define __Threading_H__

#include <AlignmentData.h>
#include <Structure.h>
#include <ThreadingInput.h>

namespace Biopool
{
/** @brief    Calculate structural scores with info derived from
*                  threading.
 * 
* @Description  
* @This 
 **/
class Threading : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	Threading(AlignmentData *ad, ThreadingInput *thread, double cThr);

	/// Copy constructor.
	Threading(const Threading &orig);

	/// Destructor.
	virtual ~Threading();


// OPERATORS:

	/// Assignment operator.
	Threading& operator = (const Threading &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const Threading &orig);

	/// Construct a new "deep copy" of this object.
	virtual Threading* newCopy();


protected:


private:

// ATTRIBUTES:

	string seq1;               ///< Target sequence.
	ThreadingInput *thread;    ///< Template threading input file.
	double cThr;               ///< Coefficient for threading.

};

} // namespace

#endif
