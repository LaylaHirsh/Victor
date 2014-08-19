// --*- C++ -*------x-----------------------------------------------------------
//
 
//
// Description:     Calculate structural scores with info derived from
//                  threading.
//
// -----------------x-----------------------------------------------------------

#include <Threading.h>

namespace Biopool
{

// CONSTRUCTORS:

Threading::Threading(AlignmentData *ad, ThreadingInput *thread, double cThr)
	: Structure(0), seq1(ad->getSequence(1)), thread(thread), cThr(cThr)
{ }


Threading::Threading(const Threading &orig) : Structure(orig)
{
	copy(orig);
}


Threading::~Threading()
{ }


// OPERATORS:

Threading&
Threading::operator = (const Threading &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
Threading::scoringStr(int i, int j)
{
	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	char aaTarget = seq1[i-1];
	int targetIndex = 0;

	for (int k = 0; k < 20; k++)
		if (aaTarget == residue_indices[k])
		{
			targetIndex = k;
			break;
		}

	return cThr * thread->score(targetIndex, (j-1));
}


// MODIFIERS:

void
Threading::copy(const Threading &orig)
{
	Structure::copy(orig);
	seq1 = orig.seq1;
	thread = orig.thread->newCopy();
	cThr = orig.cThr;
}


Threading*
Threading::newCopy()
{
	Threading *tmp = new Threading(*this);
	return tmp;
}

} // namespace
