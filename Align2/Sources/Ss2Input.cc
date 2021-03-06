// --*- C++ -*------x-----------------------------------------------------------
 
// Date:            12/2007
//
// Description:     Implement I/O objects for handling PSI-PRED files.
//
// -----------------x-----------------------------------------------------------

#include <Ss2Input.h>

namespace Biopool
{

// CONSTRUCTORS:

Ss2Input::Ss2Input()
{ }


Ss2Input::Ss2Input(istream &is)
{
	is >> *this;
}


Ss2Input::Ss2Input(const Ss2Input &orig)
{
	copy(orig);
}


Ss2Input::~Ss2Input()
{ }


// OPERATORS:

Ss2Input&
Ss2Input::operator = (const Ss2Input &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


ostream&
operator << (ostream &os, const Ss2Input &object)
{
	Ss2Input::pWriteDoubleVector(os, object.residuescores);
	return os;
}


istream&
operator >> (istream &is, Ss2Input &object)
{
	for (int i = 0; i < 8; i++) // all the first 9 words are just piped in a string variable
	{
		string tmp;
		is >> tmp;
	}

	Ss2Input::pReadDoubleVector(is, object.residuescores);
	return is;
}


// MODIFIERS:

void
Ss2Input::copy(const Ss2Input &orig)
{
	residuescores.clear();
	for (unsigned int i = 0; i < orig.residuescores.size(); i++)
	{
		vector<double> tmp;
		for (unsigned int j = 0; j < orig.residuescores[i].size(); j++)
			tmp.push_back(orig.residuescores[i][j]);
		residuescores.push_back(tmp);
	}
}


Ss2Input*
Ss2Input::newCopy()
{
	Ss2Input *tmp = new Ss2Input(*this);
	return tmp;
}


// HELPERS:

template<class T> void
Ss2Input::pWriteDoubleVector(ostream &os, vector< vector<T> > data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		for (unsigned int j = 0; j < data[j].size(); j++)
			os << data[i][j] << " ";
		os << endl;
	}
	os << endl;
}


template<class T> void
Ss2Input::pReadDoubleVector(istream &is, vector< vector<T> > &data)
{
	int tmp = -5; // setup to a funcy number; check if we have this numer

	while (is.good() and !is.eof())
	{
		int position; // 3 entries for position
		is >> position;
		if (position == tmp)
			break;
		tmp = position;

		char aa;
		is >> aa;
		char ssaa;
		is >> ssaa;

		vector<T> row;
		row.reserve(3);
		for (unsigned int j = 0; j < 3; j++) // the last 3 position for the entries in ss2 file
		{
			T tmp;
			is >> tmp;
			row.push_back(tmp);
		}
		data.push_back(row);
	}
}


template void
Ss2Input::pReadDoubleVector(istream &is, vector< vector<int> > &data);


template void
Ss2Input::pReadDoubleVector(istream &is, vector< vector<double> > &data);

} // namespace
