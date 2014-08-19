 

#ifndef __ThreadingSs2_H__
#define __ThreadingSs2_H__

#include <AlignmentData.h>
#include <Ss2Input.h>
#include <Structure.h>
#include <ThreadingInput.h>

namespace Biopool
{
/** @brief   Calculate structural scores with info derived from
*                  threading and PSI-PRED.
 * 
* @Description  
* @This 
 **/
class ThreadingSs2 : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ThreadingSs2(SubMatrix *subStr, AlignmentData *ad, ThreadingInput *thread,
		Ss2Input *psipred1, Ss2Input *psipred2, double cThr, double cSs2);

	/// Copy constructor.
	ThreadingSs2(const ThreadingSs2 &orig);

	/// Destructor.
	virtual ~ThreadingSs2();


// OPERATORS:

	/// Assignment operator.
	ThreadingSs2& operator = (const ThreadingSs2 &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ThreadingSs2 &orig);

	/// Construct a new "deep copy" of this object.
	virtual ThreadingSs2* newCopy();


protected:


private:

// ATTRIBUTES:

	string seq1;               ///< Target sequence.
	string sec1;               ///< Target secondary structure.
	string sec2;               ///< Template secondary structure.
	ThreadingInput *thread;    ///< Template threading input file.
	Ss2Input *psipred1;        ///< Target PSI-PRED input file.
	Ss2Input *psipred2;        ///< Template PSI-PRED input file.
	double cThr;               ///< Coefficient for threading.
	double cSs2;               ///< Coefficient for PSI-PRED prediction.

};

} // namespace

#endif
