// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate scores for profile to profile alignment using
//                  euclidean distance.
//
// -----------------x-----------------------------------------------------------

#include <EDistance.h>

namespace Biopool
{

// CONSTRUCTORS:

EDistance::EDistance(Profile *pro1, Profile *pro2) : ScoringFunction(),
	pro1(pro1), pro2(pro2), offset(1.00)
{ }


EDistance::EDistance(Profile *pro1, Profile *pro2, double offset)
	: ScoringFunction(), pro1(pro1), pro2(pro2), offset(offset)
{ }


EDistance::EDistance(const EDistance &orig) : ScoringFunction(orig)
{
	copy(orig);
}


EDistance::~EDistance()
{ }


// OPERATORS:

EDistance&
EDistance::operator = (const EDistance &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
EDistance::scoringSeq(int i, int j)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	double freq1 = 0.00;
	double freq2 = 0.00;
	double s = 0.00;

	for (unsigned int k = 0; k < 20; k++)
	{
		freq1 = pro1->getAminoFrequency(residue_indices[k], (i-1));
		freq2 = pro2->getAminoFrequency(residue_indices[k], (j-1));

		s += ((freq1 - freq2) * (freq1 - freq2));
	}

	return offset - sqrt(s);
}


// MODIFIERS:

void
EDistance::copy(const EDistance &orig)
{
	ScoringFunction::copy(orig);
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();
	offset = orig.offset;
}


EDistance*
EDistance::newCopy()
{
	EDistance *tmp = new EDistance(*this);
	return tmp;
}

} // namespace
