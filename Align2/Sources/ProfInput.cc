// --*- C++ -*------x-----------------------------------------------------------
//
// Description:     Implement I/O objects for handling PHD files.
//
// -----------------x-----------------------------------------------------------

#include <ProfInput.h>

namespace Biopool
{

// CONSTRUCTORS:

ProfInput::ProfInput()
{ }


ProfInput::ProfInput(istream &is)
{
	is >> *this;
}


ProfInput::ProfInput(const ProfInput &orig)
{
	copy(orig);
}


ProfInput::~ProfInput()
{ }


// OPERATORS:

ProfInput&
ProfInput::operator = (const ProfInput &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


ostream&
operator << (ostream &os, const ProfInput &object)
{
	ProfInput::pWriteString(os, object.seq, object.profSSPred,
		object.profBEPred);
	return os;
}


istream&
operator >> (istream &is, ProfInput &object)
{
	ProfInput::pReadString(is, object.seq, object.profSSPred,
		object.profBEPred);
	return is;
}


// PREDICATES:

char
ProfInput::getProfMixSSBE(char ssChar, char beChar)
{
	char correctChar;

	if (ssChar == 'H') // H case
	{
		correctChar = 'H';
		if (beChar == 'e')
			correctChar = 'I';
		return correctChar;
	}
	else
		if (ssChar == 'E') // E case
		{
			correctChar = 'E';
			if (beChar == 'e')
				correctChar = 'F';
			return correctChar;
		}
		else // C case
		{
			correctChar = 'C';
			if (beChar == 'e')
				correctChar = 'D';
			return correctChar;
		}
}


// MODIFIERS:

void
ProfInput::copy(const ProfInput &orig)
{
	seq = orig.seq;
	profSSPred = orig.profSSPred;
	profBEPred = orig.profBEPred;
}


ProfInput*
ProfInput::newCopy()
{
	ProfInput *tmp = new ProfInput(*this);
	return tmp;
}


// HELPERS:

void
ProfInput::pWriteString(ostream &os, string data1, string data2, string data3)
{
	os << "     #    AA   Pss   Pbe\n" << endl;

	for (unsigned int i = 0; i < data1.size(); i++)
	{
		os << setw(6) << i
		   << setw(6) << data1[i]
		   << setw(6) << data2[i]
		   << setw(6) << data3[i] << endl;
	}
	os << endl;
}


void
ProfInput::pReadString(istream &is, string &data1, string &data2, string &data3)
{
	string line;
	data1 = "";
	data2 = "";
	data3 = "";
	unsigned int lineType = 0;

	while ((is.good()) && (!is.eof()))
		switch(lineType)
		{
			case 0:
				line = readLine(is);
				lineType++;
				break;
			case 1:
				data1 += readLine(is);
				lineType++;
				break;
			case 2:
				data2 += readLine(is);
				lineType++;
				break;
			case 3:
				data3 += readLine(is);
				lineType = 0;
		}
}

} // namespace
