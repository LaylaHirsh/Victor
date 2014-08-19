 

#ifndef __ReverseScore_H__
#define __ReverseScore_H__

#include <Align.h>
#include <Alignment.h>
#include <StatTools.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace Biopool
{
/** @brief  
 * 
* @Description  
* @This 
 **/
class ReverseScore
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	ReverseScore(Align *a);

	/// Copy constructor.
	ReverseScore(const ReverseScore &orig);

	/// Destructor
	virtual ~ReverseScore();


// OPERATORS:

	/// Assignment operator.
	ReverseScore& operator = (const ReverseScore &orig);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const ReverseScore &orig);

	/// Calculate Z-score between ali and its reversal.
	double getZScore(double &forward, double &reverse, unsigned int n = 50);


protected:

// ATTRIBUTES:

	Align *ali;    ///< Pointer to initial Align.
	Align *inv;    ///< Pointer to inverted Align.


private:

};

} // namespace

#endif
