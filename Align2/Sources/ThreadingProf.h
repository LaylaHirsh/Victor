 

#ifndef __ThreadingProf_H__
#define __ThreadingProf_H__

#include <AlignmentData.h>
#include <ProfInput.h>
#include <Structure.h>
#include <ThreadingInput.h>

namespace Biopool
{
/** @brief    Calculate structural scores with info derived from
*                  threading and PHD.
 * 
* @Description  
* @This 
 **/
class ThreadingProf : public Structure
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ThreadingProf(SubMatrix *subStr, AlignmentData *ad, ThreadingInput *thread,
		ProfInput *phd1, ProfInput *phd2, double cThr, double cPrf);

	/// Copy constructor.
	ThreadingProf(const ThreadingProf &orig);

	/// Destructor.
	virtual ~ThreadingProf();


// OPERATORS:

	/// Assignment operator.
	ThreadingProf& operator = (const ThreadingProf &orig);


// PREDICATES:

	/// Calculate structural scores to create matrix values.
	virtual double scoringStr(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ThreadingProf &orig);

	/// Construct a new "deep copy" of this object.
	virtual ThreadingProf* newCopy();


protected:


private:

// ATTRIBUTES:

	string seq1;               ///< Target sequence.
	ThreadingInput *thread;    ///< Template threading input file.
	ProfInput *phd1;           ///< Target PHD input file.
	ProfInput *phd2;           ///< Template PHD input file.
	double cThr;               ///< Coefficient for threading.
	double cPrf;               ///< Coefficient for PROF prediction.

};

} // namespace

#endif
