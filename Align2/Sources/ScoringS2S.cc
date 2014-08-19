// --*- C++ -*------x-----------------------------------------------------------
//
// Description:     Calculate scores for sequence to sequence alignment.
//
// -----------------x-----------------------------------------------------------

#include <ScoringS2S.h>

namespace Biopool
{

// CONSTRUCTORS:

ScoringS2S::ScoringS2S(SubMatrix *sub, AlignmentData *ad, Structure *str,
	double cSeq) : ScoringScheme(sub, ad, str), seq1(ad->getSequence(1)),
		seq2(ad->getSequence(2)), cSeq(cSeq)
{ }


ScoringS2S::ScoringS2S(const ScoringS2S &orig) : ScoringScheme(orig)
{
	copy(orig);
}


ScoringS2S::~ScoringS2S()
{ }


// OPERATORS:

ScoringS2S&
ScoringS2S::operator = (const ScoringS2S &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
ScoringS2S::scoring(int i, int j)
{
	double s = cSeq * sub->score[seq1[i-1]][seq2[j-1]];

	if (str != 0)
		s += str->scoringStr(i, j);

	return s;
}


// MODIFIERS:

void
ScoringS2S::copy(const ScoringS2S &orig)
{
	ScoringScheme::copy(orig);
	seq1 = orig.seq1;
	seq2 = orig.seq2;
	cSeq = orig.cSeq;
}


ScoringS2S*
ScoringS2S::newCopy()
{
	ScoringS2S *tmp = new ScoringS2S(*this);
	return tmp;
}


void
ScoringS2S::reverse()
{
	ScoringScheme::reverse();

	string tmp = "";
	for (unsigned int i = seq2.length(); i > 0; i--)
		tmp.push_back(seq2[i-1]);
	seq2 = tmp;
}

} // namespace
