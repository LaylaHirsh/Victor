// --*- C++ -*------x-----------------------------------------------------------
//
// Description:     Calculate scores for profile to profile alignment.
//
// -----------------x-----------------------------------------------------------

#include <ScoringP2P.h>

namespace Biopool
{

// CONSTRUCTORS:

ScoringP2P::ScoringP2P(SubMatrix *sub, AlignmentData *ad, Structure *str,
	Profile *pro1, Profile *pro2, ScoringFunction *fun, double cSeq)
		: ScoringScheme(sub, ad, str), seq1(ad->getSequence(1)),
			seq2(ad->getSequence(2)), pro1(pro1), pro2(pro2), fun(fun),
				cSeq(cSeq)
{ }


ScoringP2P::ScoringP2P(const ScoringP2P &orig) : ScoringScheme(orig)
{
	copy(orig);
}


ScoringP2P::~ScoringP2P()
{ }


// OPERATORS:

ScoringP2P&
ScoringP2P::operator = (const ScoringP2P &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
ScoringP2P::scoring(int i, int j)
{
	//cout<<"scoring A\n";
	double s = cSeq * fun->scoringSeq(i, j);
	//cout<<"scoring B\n";
	if (str != 0)
		s += str->scoringStr(i, j);
	//cout<<"scoring C\n";
	return s;
}


// MODIFIERS:

void
ScoringP2P::copy(const ScoringP2P &orig)
{
	ScoringScheme::copy(orig);
	seq1 = orig.seq1;
	seq2 = orig.seq2;
	pro1 = orig.pro1->newCopy();
	pro2 = orig.pro2->newCopy();
	fun = orig.fun->newCopy();
	cSeq = orig.cSeq;
}


ScoringP2P*
ScoringP2P::newCopy()
{
	ScoringP2P *tmp = new ScoringP2P(*this);
	return tmp;
}


void
ScoringP2P::reverse()
{
	ScoringScheme::reverse();

	string tmp = "";
	for (unsigned int i = seq2.length(); i > 0; i--)
		tmp.push_back(seq2[i-1]);
	seq2 = tmp;

	pro2->reverse();
}

} // namespace
