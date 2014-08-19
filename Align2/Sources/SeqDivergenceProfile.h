 
#ifndef __SeqDivergenceProfile_H__
#define __SeqDivergenceProfile_H__

#include <AGPFunction.h>
#include <NWAlign.h>
#include <Profile.h>
#include <ScoringS2S.h>
#include <SequenceData.h>
#include <SubMatrix.h>
#include <math.h>

namespace Biopool
{/** @brief    Calculate a frequency profile or PSSM using SeqDivergence
*                  weighting scheme.
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
*
* @This 
 **/

class SeqDivergenceProfile : public Profile
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	SeqDivergenceProfile();

	/// Copy constructor.
	SeqDivergenceProfile(const SeqDivergenceProfile &orig);

	/// Destructor.
	virtual ~SeqDivergenceProfile();


// OPERATORS:

	/// Assignment operator.
	SeqDivergenceProfile& operator = (const SeqDivergenceProfile &orig);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const SeqDivergenceProfile &orig);

	/// Construct a new "deep copy" of this object.
	virtual SeqDivergenceProfile* newCopy();


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

	vector<double> aliWeight;    ///< Alignment weights.

};

} // namespace

#endif
