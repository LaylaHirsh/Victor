// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate scores for profile to profile alignment using
//                  dot product method. Some explanations can be found in:
//
// 	                Mittelman D., Sadreyev R., Grishin N.
//                  Probabilistic scoring measures for profile-profile
//                  comparison yield more accurate short seed alignments.
//                  Bioinformatics. 2003 Aug 12;19(12):1531-9.
//                  PMID: 12912834 [PubMed - indexed for MEDLINE]
//
//                  Marti-Renom MA., Madhusudhan MS., Sali A.
//                  Alignment of protein sequences by their profiles.
//                  Protein Sci. 2004 Apr;13(4):1071-87.
//                  PMID: 15044736 [PubMed - indexed for MEDLINE]
//
// -----------------x-----------------------------------------------------------

#include <DotPFreq.h>

namespace Biopool
{

// CONSTRUCTORS:

DotPFreq::DotPFreq(Profile *pro1, Profile *pro2) : ScoringFunction(),
	pro1(pro1), pro2(pro2)
{ }


DotPFreq::DotPFreq(const DotPFreq &orig) : ScoringFunction(orig)
{
	copy(orig);
}


DotPFreq::~DotPFreq()
{ }


// OPERATORS:

DotPFreq&
DotPFreq::operator = (const DotPFreq &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
DotPFreq::scoringSeq(int i, int j)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	double freq1 = 0.00;
	double freq2 = 0.00;
	double s = 0.00;

	for (unsigned int k = 0; k < 20; k++)
	{
		freq1 = pro1->getAminoFrequency(residue_indices[k], (i-1));
		freq2 = pro2->getAminoFrequency(residue_indices[k], (j-1));

		s += (freq1 * freq2);
	}

	return s;
}


// MODIFIERS:

void
DotPFreq::copy(const DotPFreq &orig)
{
	ScoringFunction::copy(orig);
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();
}


DotPFreq*
DotPFreq::newCopy()
{
	DotPFreq *tmp = new DotPFreq(*this);
	return tmp;
}

} // namespace
