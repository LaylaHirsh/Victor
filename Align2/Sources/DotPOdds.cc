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

#include <DotPOdds.h>

namespace Biopool
{

// CONSTRUCTORS:

DotPOdds::DotPOdds(Profile *pro1, Profile *pro2) : ScoringFunction(),
	pro1(pro1), pro2(pro2)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	for (unsigned int k = 0; k < 20; k++)
	{
		p1[k] = 0.00001;
		for (unsigned int n = 0; n < pro1->getSequenceLength(); n++)
			p1[k] += (pro1->getAminoFrequency(residue_indices[k], n) / pro1->getSequenceLength());
	}

	for (unsigned int k = 0; k < 20; k++)
	{
		p2[k] = 0.00001;
		for (unsigned int n = 0; n < pro2->getSequenceLength(); n++)
			p2[k] += (pro2->getAminoFrequency(residue_indices[k], n) / pro2->getSequenceLength());
	}
}


DotPOdds::DotPOdds(const DotPOdds &orig) : ScoringFunction(orig)
{
	copy(orig);
}


DotPOdds::~DotPOdds()
{ }


// OPERATORS:

DotPOdds&
DotPOdds::operator = (const DotPOdds &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
DotPOdds::scoringSeq(int i, int j)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	double freq1 = 0.00;
	double freq2 = 0.00;
	double odds1 = 0.00;
	double odds2 = 0.00;
	double s = 0.00;

	for (unsigned int k = 0; k < 20; k++)
	{
		freq1 = pro1->getAminoFrequency(residue_indices[k], (i-1)) + 0.00001;
		freq2 = pro2->getAminoFrequency(residue_indices[k], (j-1)) + 0.00001;

		odds1 = log(freq1/p1[k]);
		odds2 = log(freq2/p2[k]);

		s += (odds1 * odds2);
	}

	return s;
}


// MODIFIERS:

void
DotPOdds::copy(const DotPOdds &orig)
{
	ScoringFunction::copy(orig);
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();

	double p1[20];
	for (unsigned int k = 0; k < 20; k++)
		p1[k] = orig.p1[k];

	double p2[20];
	for (unsigned int k = 0; k < 20; k++)
		p2[k] = orig.p2[k];
}


DotPOdds*
DotPOdds::newCopy()
{
	DotPOdds *tmp = new DotPOdds(*this);
	return tmp;
}

} // namespace
