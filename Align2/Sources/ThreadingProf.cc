// --*- C++ -*------x-----------------------------------------------------------
 
// Description:     Calculate structural scores with info derived from
//                  threading and PHD.
//
// -----------------x-----------------------------------------------------------

#include <ThreadingProf.h>

namespace Biopool
{

// CONSTRUCTORS:

ThreadingProf::ThreadingProf(SubMatrix *subStr, AlignmentData *ad,
	ThreadingInput *thread, ProfInput *phd1, ProfInput *phd2, double cThr,
		double cPrf) : Structure(subStr), seq1(ad->getSequence(1)),
			thread(thread), phd1(phd1), phd2(phd2), cThr(cThr), cPrf(cPrf)
{ }


ThreadingProf::ThreadingProf(const ThreadingProf &orig) : Structure(orig)
{
	copy(orig);
}


ThreadingProf::~ThreadingProf()
{ }


// OPERATORS:

ThreadingProf&
ThreadingProf::operator = (const ThreadingProf &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
ThreadingProf::scoringStr(int i, int j)
{
	//
	// THREADING
	//

	const string residue_indices = "ARNDCQEGHILKMFPSTWYV";

	char aaTarget = seq1[i-1];
	int targetIndex = 0;

	for (int k = 0; k <= 19; k++)
		if (aaTarget == residue_indices[k])
		{
			targetIndex = k;
			break;
		}

	double s1 = thread->score(targetIndex, (j-1));


	//
	// PROF
	//

	// PROFsec prediction
	char ssTarget = phd1->getProfSS(i-1);
	char ssTemplate = phd2->getProfSS(j-1);

	// PROFacc prediction
	char beTarget = phd1->getProfBE(i-1);
	char beTemplate = phd2->getProfBE(j-1);

	char mixTarget = phd1->getProfMixSSBE(ssTarget, beTarget);
	char mixTemplate = phd2->getProfMixSSBE(ssTemplate, beTemplate);

	double s2 = subStr->score[mixTarget][mixTemplate];


	return cThr * s1 + cPrf * s2;
}


// MODIFIERS:

void
ThreadingProf::copy(const ThreadingProf &orig)
{
	Structure::copy(orig);
	seq1 = orig.seq1;
	thread = orig.thread->newCopy();
	phd1 = orig.phd1->newCopy();
	phd2 = orig.phd2->newCopy();
	cThr = orig.cThr;
	cPrf = orig.cPrf;
}


ThreadingProf*
ThreadingProf::newCopy()
{
	ThreadingProf *tmp = new ThreadingProf(*this);
	return tmp;
}

} // namespace
