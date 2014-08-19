// --*- C++ -*------x-----------------------------------------------------------
//
//
// Description:     Calculate structural scores with info derived from PHD.
//
// -----------------x-----------------------------------------------------------

#include <Prof.h>

namespace Biopool
{

// CONSTRUCTORS:

Prof::Prof(SubMatrix *subStr, ProfInput *phd1, ProfInput *phd2, double cPrf)
	: Structure(subStr), phd1(phd1), phd2(phd2), cPrf(cPrf)
{ }


Prof::Prof(const Prof &orig) : Structure(orig)
{
	copy(orig);
}


Prof::~Prof()
{ }


// OPERATORS:

Prof&
Prof::operator = (const Prof &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

double
Prof::scoringStr(int i, int j)
{
	// ProfSEC prediction
	char ssTarget = phd1->getProfSS(i-1);
	char ssTemplate = phd2->getProfSS(j-1);

	// ProfASA prediction
	char beTarget = phd1->getProfBE(i-1);
	char beTemplate = phd2->getProfBE(j-1);

	char mixTarget = phd1->getProfMixSSBE(ssTarget, beTarget);
	char mixTemplate = phd2->getProfMixSSBE(ssTemplate, beTemplate);

	return cPrf * subStr->score[mixTarget][mixTemplate];
}


// MODIFIERS:

void
Prof::copy(const Prof &orig)
{
	Structure::copy(orig);
	phd1 = orig.phd1->newCopy();
	phd2 = orig.phd2->newCopy();
	cPrf = orig.cPrf;
}


Prof*
Prof::newCopy()
{
	Prof *tmp = new Prof(*this);
	return tmp;
}

} // namespace
