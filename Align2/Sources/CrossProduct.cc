// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate scores for profile to profile alignment using
//                  sum of pairs method. Some explanations can be found in:
//
//                  Guoli Wang, Roland L. Dunbrack jr.
//                  Scoring profile-to-profile sequence alignments.
//                  Institute for Cancer Research, Fox Chase Cancer Center,
//                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
//
// -----------------x-----------------------------------------------------------

              /*************************************************
               * sub is BLOSUM62 matrix standard log-odds form *
               *************************************************/

#include <CrossProduct.h>

namespace Biopool
{

// CONSTRUCTORS:

CrossProduct::CrossProduct(SubMatrix *sub, Profile *pro1, Profile *pro2)
	: ScoringFunction(), sub(sub), pro1(pro1), pro2(pro2)
{ }


CrossProduct::CrossProduct(const CrossProduct &orig) : ScoringFunction(orig)
{
	copy(orig);
}


CrossProduct::~CrossProduct()
{ }


// OPERATORS:

CrossProduct&
CrossProduct::operator = (const CrossProduct &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
CrossProduct::scoringSeq(int i, int j)
{
	double s = 0.00;

	for (AminoAcidCode amino1 = ALA; amino1 <= TYR; amino1++)
	{
		double tmp = 0.00;

		for (AminoAcidCode amino2 = ALA; amino2 <= TYR; amino2++)
			tmp += sub->score[aminoAcidOneLetterTranslator(amino1)]
				[aminoAcidOneLetterTranslator(amino2)] *
					pro2->getAminoFrequencyFromCode(amino2, (j-1));

		s += (pro1->getAminoFrequencyFromCode(amino1, (i-1)) * tmp);
	}

	return s;
}


// MODIFIERS:

void
CrossProduct::copy(const CrossProduct &orig)
{
	ScoringFunction::copy(orig);
	sub = orig.sub->newCopy();
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();
}


CrossProduct*
CrossProduct::newCopy()
{
	CrossProduct *tmp = new CrossProduct(*this);
	return tmp;
}

} // namespace
