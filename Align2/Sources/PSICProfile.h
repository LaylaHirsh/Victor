 

#ifndef __PSICProfile_H__
#define __PSICProfile_H__

#include <Profile.h>

namespace Biopool
{
/** @brief   Calculate a frequency profile or PSSM using PSIC weighting
*                  scheme.
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class PSICProfile : public Profile
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	PSICProfile();

	/// Copy constructor.
	PSICProfile(const PSICProfile &orig);

	/// Destructor.
	virtual ~PSICProfile();


// OPERATORS:

	/// Assignment operator.
	PSICProfile& operator = (const PSICProfile &orig);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const PSICProfile &orig);

	/// Construct a new "deep copy" of this object.
	virtual PSICProfile* newCopy();


// HELPERS:

	/// Calculate alignment weights.
	virtual void pCalculateWeight(Alignment &ali);

	/// Calculate the raw (ie. unnormalized) aminoacids frequencies for position i.
	virtual void pCalculateRawFrequency(vector<double> &freq, double &gapFreq,
		Alignment &ali, unsigned int i);

	/// Construct data from alignment.
	virtual void pConstructData(Alignment &ali);


protected:


private:

	vector< vector<double> > aliWeight;    ///< Alignment weights.

};

} // namespace

#endif
