// --*- C++ -*------x-----------------------------------------------------------
//
 
// Description:     Base class for scoring schemes.
//
// -----------------x-----------------------------------------------------------

#include <ScoringScheme.h>

namespace Biopool
{

// CONSTRUCTORS:

ScoringScheme::ScoringScheme(SubMatrix *sub, AlignmentData *ad, Structure *str)
	: sub(sub), ad(ad), str(str)
{
	if ((!checkSequence(ad->getSequence(1))) ||
		(!checkSequence(ad->getSequence(2))))
	{
		cout << "Illegal sequence:\n"
		     << ad->getSequence(1) << "\n"
		     << ad->getSequence(2) << "\n"
		     << sub->getResidues() << endl;
		ERROR("Error checking sequence.", exception);
	}
}


ScoringScheme::ScoringScheme(const ScoringScheme &orig)
{
	copy(orig);
}


ScoringScheme::~ScoringScheme()
{ }


// OPERATORS:

ScoringScheme&
ScoringScheme::operator = (const ScoringScheme &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

bool
ScoringScheme::checkSequence(const string &s) const
{
	string residues = sub->getResidues();
	bool found = false;

	for (unsigned int i = 0; i < s.size(); i++)
	{
		found = false;

		for (unsigned int j = 0; j < residues.size(); j++)
			if (residues[j] == s[i])
			{
				found = true;
				break;
			}

		if (!found)
		{
			cout << "checkSequence unsuccessful!" << endl;
			return false;
		}
	}

	return true;
}


// MODIFIERS:

void
ScoringScheme::copy(const ScoringScheme &orig)
{
	sub = orig.sub->newCopy();
	ad = orig.ad->newCopy();
	str = orig.str->newCopy();
}


void
ScoringScheme::reverse()
{
	if (str != 0)
		str->reverse();
}

} // namespace
